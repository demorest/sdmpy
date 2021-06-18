from __future__ import print_function, division, absolute_import, unicode_literals # not casa compatible
from builtins import bytes, dict, object, range, map, input#, str # not casa compatible
from future.utils import itervalues, viewitems, iteritems, listvalues, listitems
from io import open

import os
import sys
import re
import mmap
import math
import numpy
from copy import deepcopy
from lxml import etree, objectify

try:
    from progressbar import ProgressBar
except ImportError:
    ProgressBar = None

from .mime import MIMEPart, MIMEHeader

import logging
logger = logging.getLogger(__name__)


# TODO find a better way to get the namespace automatically?
_ns = '{http://Alma/XASDM/sdmbin}'


def basename_noext(path):
    return os.path.basename(os.path.splitext(path)[0])


def _stripns(tag):
    return re.sub('{.+}', '', tag)


def ant2bl(i, j=None):
    """Returns baseline index for given antenna pair.  Will accept
    two args, or a list/tuple/etc.  Uses 0-based indexing"""
    if j is None:
        (a1, a2) = sorted(i[:2])
    else:
        (a1, a2) = sorted((i, j))
    # could raise error if a2==a1, either are negative, etc
    return (a2*(a2-1))//2 + a1


def bl2ant(i):
    """Returns antenna pair for given baseline index.  All are 0-based."""
    a2 = int(0.5*(1.0+math.sqrt(1.0+8.0*i)))
    a1 = i - a2*(a2-1)//2
    return a1, a2


# This class provides access to EVLA/ALMA Binary Data Format (BDF)
# files.  The approach used here is based loosely on 'bdfparse.py'
# originally by Peter Williams.  So far this is somewhat targeted
# to EVLA/WIDAR data files, it may not handle all the ALMA variants.

_mmap_default = 'auto'
_mmap_limit = 8<<30

class BDF(object):
    """
    Class representing a single BDF file.  For example:

        b = bdf.BDF('uid____evla_bdf_1433189755525')

    Individual integration data is returned as BDFIntegration objects via
    either b.get_integration(idx) or b[idx].  Other useful methods include:

        b.basebands      # list of baseband ids
        b.spws           # list of spectral windows per baseband
        b.numAntenna     # number of antennas
        b.numBaseline    # number of baselines
        b.numIntegration # number of integrations in file
        b.sdmDataHeader  # lxml objectify version of full header

    The constructor takes a kwarg, use_mmap to specify whether the data
    are read using mmap versus a simple read().  The former can handle
    arbitrarily large files easily but seems to have poor performance on
    some filesystems (lustre).  If use_mmap=False is specified you must
    ensure there is enough free memory to hold the entire BDF file.
    Allowed values for use_mmap are:
    
        True:       always use mmap
        False:      never use mmap
        'auto':     use mmap for files larger than sdmpy.bdf._mmap_limit
                    (default 8GB) 
        'default':  Apply setting from sdmpy.bdf._mmap_default

    """

    def __init__(self, fname, use_mmap='default'):
        self.fname = fname
        try:
            self.fp = open(fname, 'rb')
        except IOError:
            self.fp = None
        else:
            if use_mmap=='default':
                use_mmap = _mmap_default
            if use_mmap=='auto':
                fsize = self.fp.seek(0,os.SEEK_END)
                self.fp.seek(0)
                use_mmap = (fsize > _mmap_limit)
            if use_mmap:
                self.mmdata = mmap.mmap(self.fp.fileno(), 0, mmap.MAP_PRIVATE,
                        mmap.PROT_READ)
            else:
                self.mmdata = self.fp.read()
            self.read_mime()
            self.parse_spws()

    # Size in bytes of each data element. In principle, the crossData
    # type needs to be read from the headers while all others have
    # pre-set values...
    bin_dtype_size = {
            'flags':           4,  # INT32
            'actualTimes':     8,  # INT64
            'actualDurations': 8,  # INT64
            'zeroLags':        4,  # FLOAT32
            'autoData':        4,  # FLOAT32
            'crossData':       4,  # FLOAT32 but should be determined from file
            }

    # Basic data type for each array.  Note, does not necessarily
    # match up with the sizes above due to how BDF defines the 'size'
    # attribute..
    bin_dtype = {
            'flags':           numpy.int32,
            'actualTimes':     numpy.int64,
            'actualDurations': numpy.int64,
            'zeroLags':        numpy.float32,
            'autoData':        numpy.float32,
            'crossData':       numpy.complex64,
            }

    @property
    def exists(self):
        return self.fp is not None

    def read_mime(self, full_read=False):
        if self.fp:
            self.fp.seek(0, 0)  # Go back to start
            if not self.fp.readline().decode('utf-8').startswith('MIME-Version:'):
                raise RuntimeError('Invalid BDF: missing MIME-Version')

            # First we need to read and parse only the main XML header in order
            # to get sizes of the binary parts.  Note, the info stored in
            # self.bin_size is in bytes, rather than the weird BDF units.
            mime_hdr = MIMEPart(self.fp).hdr
            self.top_mime_bound = mime_hdr.boundary
            sdmDataMime = MIMEPart(self.fp, boundary=self.top_mime_bound)
            if sdmDataMime.loc != 'sdmDataHeader.xml':
                raise RuntimeError('Invalid BDF: missing sdmDataHeader.xml')
            self.sdmDataHeader = objectify.fromstring(bytes(sdmDataMime.body, 'utf-8'))
            self.bin_size = {}
            self.bin_axes = {}
            for e in self.sdmDataHeader.iter():
                if 'size' in list(e.attrib.keys()) and 'axes' in list(e.attrib.keys()):
                    binname = _stripns(e.tag)
                    self.bin_size[binname] = int(e.attrib['size']) \
                        * self.bin_dtype_size[binname]
                    self.bin_axes[binname] = e.attrib['axes'].split()

            # For EVLA, we can read the first integration, note the file offset
            # in order to determine the size, then seek to each integration as
            # requested, rather than parsing the whole file here.
            if 'EVLA' in mime_hdr['Content-Description'][0] and not full_read:
                self.offset_ints = self.fp.tell()  # Offset in file to first integ
                self.mime_ints = [MIMEPart(self.fp,
                                           boundary=self.top_mime_bound,
                                           binary_size=self.bin_size,
                                           recurse=True), ]
            # Compute size of each integration section:
                self.size_ints = self.fp.tell() - self.offset_ints
                numints = int((os.path.getsize(self.fname)-self.offset_ints)//self.size_ints)
                self.mime_ints += [None, ]*(numints-1)

            # This is the more general way to do it that does not assume
            # each integration (including XML and MIME headers) has the
            # same size in the file.  In this case, go back to the beginning
            # and parse the whole MIME structure to map it out.
            else:
                self.fp.seek(0, 0)  # reset to start again
                full_mime = MIMEPart(self.fp,
                                     recurse=True, binary_size=self.bin_size)
                self.mime_ints = full_mime.body[1:]
        else:
            logger.warn('No BDF file found at {0}'.format(self.fname))


    def _raw(self, idx):
        if self.fp:
            if self.mime_ints[idx] is not None:
                    return self.mime_ints[idx]

            # Need to read this one
            self.fp.seek(self.offset_ints + idx*self.size_ints, 0)
            #self.mime_ints[idx] = self.read_mime_part(boundary=self.top_mime_bound,recurse=True)
            self.mime_ints[idx] = MIMEPart(self.fp,
                                           boundary=self.top_mime_bound,
                                           binary_size=self.bin_size,
                                           recurse=True)
            return self.mime_ints[idx]
        else:
            logger.warn('No BDF file found at {0}'.format(self.fname))
            return None

    @property
    def projectPath(self):
        return self.sdmDataHeader.attrib['projectPath']

    @property
    def numIntegration(self):
        return len(self.mime_ints)

    @property
    def numAntenna(self):
        return int(self.sdmDataHeader.numAntenna)

    @property
    def numBaseline(self):
        return (self.numAntenna*(self.numAntenna-1))//2

    @property
    def startTime(self):
        return float(self.sdmDataHeader.startTime)/86400.0e9

    def parse_spws(self):
        self.basebands = []
        self.spws = []
        # Offsets into the cross and data arrays where this spw
        # can be found.  This info could be reconstructed later
        # but seems convenient to do it here.
        cross_offset = 0
        auto_offset = 0
        cross_meta_offset = 0
        auto_meta_offset = 0
        self.spws = []
        for bb in self.sdmDataHeader.dataStruct.baseband:
            bbname = bb.attrib['name']
            self.basebands.append(bbname)
            for spw_elem in bb.spectralWindow:
                # Build a list of spectral windows for each baseband
                spw = BDFSpectralWindow(spw_elem, cross_offset, auto_offset,
                        cross_meta_offset=cross_meta_offset,
                        auto_meta_offset=auto_meta_offset)
                cross_offset += spw.dsize('cross')
                auto_offset += spw.dsize('auto')
                cross_meta_offset += spw.msize('cross')
                auto_meta_offset += spw.msize('auto')
                self.spws.append(spw)

    def get_integration(self, idx):
        return BDFIntegration(self, idx)

    def __getitem__(self, idx):
        return self.get_integration(idx)

    def zerofraction(self, spwidx='all', type='cross'):
        """
        Return zero fraction for the entire BDF.  This is done by loading
        each integration's data so may take a while.
        """

        tot = 0.0
        for i in self:
            tot += i.zerofraction(spwidx, type)
        return tot / self.numIntegration

    def get_data(self, spwidx='all', type='cross', scrunch=False,
                 fscrunch=False, frange=None, trange=None, bar=False):
        """Returns an array containing all integrations for the specified
        spw and data type.  Takes a number of options:

          trange: tuple giving range of integrations to return (default=all).
          scrunch: if True, all requested integrations will be time-averaged.
          fscrunch: if True, data will be averaged over the channel axis.
          frange: range of channels to average if using fscrunch.
          bar: if True and the progressbar package is available, display
            a progress bar as data are loaded.

        If spwidx=='all' and no averaging was requested, then the dimensions
        of the returned array are: (time, baseline/antenna, spw, bin, channel,
        polarization).  If a single spw is selected, the spw axis is omitted.
        If time and/or freq averaging is selected, the time and/or channel
        axes are omitted.
        """
        chidx = -2  # index of spectral channels
        # Figure out ranges:
        if trange is None:
            i0 = 0
            i1 = self.numIntegration
        else:
            i0 = trange[0]
            i1 = trange[1]
        # Read first integration to get shapes, etc
        nsubout = i1 - i0
        subdat = self.get_integration(i0).get_data(spwidx, type)
        if scrunch:
            dshape = subdat.shape
        else:
            dshape = (nsubout,) + subdat.shape
        if fscrunch:
            dshape = dshape[:chidx] + dshape[chidx+1:]
        result = numpy.zeros(dshape, dtype=subdat.dtype)
        if bar and ProgressBar is not None:
            b = ProgressBar()
        else:
            b = lambda x: x
        for i in b(list(range(i0, i1))):
            if fscrunch:
                if frange is None:
                    dat = self.get_integration(i).get_data(spwidx, type).mean(chidx)
                else:
                    dat = self.get_integration(i).get_data(spwidx,
                                                           type).take(list(range(*frange)),
                                                                      axis=chidx).mean(chidx)
            else:
                dat = self.get_integration(i).get_data(spwidx, type)
            if scrunch:
                result += dat
            else:
                result[i-i0] = dat
        if scrunch:
            result /= float(nsubout)
        return result

    def get_meta(self, component, spwidx='all', corr='cross', trange=None):
        if trange is None:
            i0 = 0
            i1 = self.numIntegration
        else:
            i0 = trange[0]
            i1 = trange[1]
        nsubout = i1 - i0
        subdat = self.get_integration(i0).get_meta(component, spwidx, corr)
        dshape = (nsubout,) + subdat.shape
        result = numpy.zeros(dshape, dtype=subdat.dtype)
        for i in range(i0, i1):
            dat = self.get_integration(i).get_meta(component, spwidx, corr)
            result[i-i0] = dat
        return result

class BDFSpectralWindow(object):
    """Class that represents spectral window information present in BDF files,
    including storing appropriate offsets into the main data array.  Should be
    initialized from the spectralWindow XML element from the main BDF
    header.

    Alternatively, a BDFSpectralWindow can be generated directly without
    reference to an existing XML element.  In this case, the following
    arguments should be filled in appropriately (names are same as in XML):
        numBin
        numSpectralPoint
        sw
        swbb
    And the npol argument should be set to either 2 or 4; other values are
    not currently handled.

    An XML Element can be generated (eg for output) using the to_xml()
    method.
    """

    # Example spw element:
    # <spectralWindow sw="1" swbb="AC_8BIT" sdPolProducts="RR RL LL" crossPolProducts="RR RL LR LL" numSpectralPoint="64" numBin="40" scaleFactor="1.000000" sideband="NOSB"/>

    def __init__(self, spw_elem, cross_offset=None, auto_offset=None,
                 numBin=None, numSpectralPoint=None, sw=None, swbb=None,
                 npol=None, cross_meta_offset=None, auto_meta_offset=None):
        if spw_elem is not None:
            self._attrib = spw_elem.attrib
        else:
            # Fill _attrib based on extra input:
            self._attrib = {}
            self._attrib['numBin'] = '%d' % numBin
            self._attrib['numSpectralPoint'] = '%d' % numSpectralPoint
            self._attrib['sw'] = '%d' % sw
            self._attrib['swbb'] = str(swbb)
            if npol == 4:
                self._attrib['sdPolProducts'] = 'RR RL LL'
                self._attrib['crossPolProducts'] = 'RR RL LR LL'
            elif npol == 2:
                self._attrib['sdPolProducts'] = 'RR LL'
                self._attrib['crossPolProducts'] = 'RR LL'
            else:
                raise RuntimeError("Don't know how to handle npol=%d" % npol)
            # Boilerplate for VLA?
            self._attrib['scaleFactor'] = '1.000000'
            self._attrib['sideband'] = 'NOSB'

        self.cross_offset = cross_offset
        self.auto_offset = auto_offset
        self.cross_meta_offset = cross_meta_offset
        self.auto_meta_offset = auto_meta_offset


    def to_xml(self):
        # Property?
        # Note this returns a standalone xml Element, _not_ a reference to
        # the original document structure that this was derived from.
        result = etree.Element('spectralWindow')
        for k, v in list(self._attrib.items()):
            result.attrib[k] = v
        return result

    @property
    def numBin(self):
        return int(self._attrib['numBin'])

    @property
    def numSpectralPoint(self):
        return int(self._attrib['numSpectralPoint'])

    # These two are not listed in the BDF spec but appear in EVLA data:
    @property
    def sw(self):
        return int(self._attrib['sw'])

    @property
    def swbb(self):
        return self._attrib['swbb']

    # Should uniquely identify a sw?
    @property
    def name(self):
        return self.swbb + '-' + str(self.sw)

    def pols(self, type):
        """Return number of polarization array elements for the given data
        type (cross or auto)."""
        try:
            if type[0].lower() == 'c':
                return self._attrib['crossPolProducts'].split()
            elif type[0].lower() == 'a':
                return self._attrib['sdPolProducts'].split()
        except KeyError:
            return 0

    def npol(self, type):
        """Return number of polarization array elements for the given data
        type (cross or auto)."""
        pols = self.pols(type=type)
        if type[0].lower() == 'c':
            return len(pols)
        elif type[0].lower() == 'a':
            return 4 if len(pols) == 3 else len(pols)  # 3==4 in BDF math! :)

    def dshape(self, type):
        """Return shape tuple of data array for this spectral window,
        in number of data elements (real for auto, complex for cross).
        """

        return (self.numBin, self.numSpectralPoint, self.npol(type))

    def dsize(self, type):
        """Return size of data array for this spectral window, in number of
        data elements (real for auto, complex for cross)."""
        return numpy.product(self.dshape(type))

    def mshape(self, corr):
        """Return shape tuple of metadata array for this spectral window
        and correlation type (cross or auto)."""
        return (self.numBin, len(self.pols(corr)))

    def msize(self, corr):
        """Return size of metadata array for this spectral window and
        correlation type (cross or auto)."""
        return numpy.product(self.mshape(corr))

    @staticmethod
    def dims_match(spwlist, type):
        """Given a list of BDFSpectralWindow objects, return true if all
        of them have consistent array dimensions."""
        if len(spwlist) == 1:
            return True
        for spw in spwlist[1:]:
            if spwlist[0].dshape(type) != spw.dshape(type):
                return False
        return True


class BDFIntegration(object):
    """
    Describes and holds data for a single intgration within a BDF file.
    This should be derived from an existing BDF object using
    get_integration() or indexing, ie:

        b = bdf.BDF('some_file')

        # Get the 5th integration, these two are equivalent:
        i = b.get_integration(5)
        i = b[5]

        # Read the cross-corr data array for spectral window 0
        dat = i.get_data(0)

    Other potentially useful info:

        i.basebands            # list of baseband IDs
        i.spws                 # dict of spws per baseband
        i.numAntenna           # obvious
        i.numBaseline          # "
        i.sdmDataSubsetHeader  # lxml objectify version of full sub-header

    """

    def __init__(self, bdf, idx):
        # Get the main header
        self.sdmDataSubsetHeader = objectify.fromstring(
                bytes(bdf._raw(idx).body[0].body, 'utf-8'))
        # Copy some info from the BDF headers
        self.basebands = bdf.basebands
        self.spws = bdf.spws
        self.bin_axes = bdf.bin_axes
        self.numAntenna = bdf.numAntenna
        self.numBaseline = bdf.numBaseline
        # Get the binary data
        self.data = {}
        for m in bdf._raw(idx).body[1:]:
            btype = basename_noext(m.loc)
            bsize = bdf.bin_size[btype]  # size of the binary blob in bytes
            baxes = self.bin_axes[btype]
            # Determine outer step size for the array, either baselines
            # antennas or baseline+antenna.  We can't apply the other
            # dimensions here because the number of elements can vary
            # per spw.
            if baxes[0] == 'BAL' and baxes[1] == 'ANT':
                # In this case, cross and auto may have different step sizes
                # so nothing we can do here.
                shape = (-1,)
            elif baxes[0] == 'BAL':
                shape = (self.numBaseline, -1)
            elif baxes[0] == 'ANT':
                shape = (self.numAntenna, -1)
            else:
                shape = (-1,)  # Don't know what to do, just leave flat array
            self.data[btype] = numpy.frombuffer(bdf.mmdata[m.body:m.body+bsize],
                                                dtype=bdf.bin_dtype[btype]).reshape(shape)

    @property
    def projectPath(self):
        return self.sdmDataSubsetHeader.attrib['projectPath']

    @property
    def time(self):
        return float(self.sdmDataSubsetHeader.schedulePeriodTime.time)/86400.0e9

    @property
    def interval(self):
        return float(self.sdmDataSubsetHeader.schedulePeriodTime.interval)*1e-9

    def get_data(self, spwidx='all', type='cross'):
        """
        Return the data array for the given subset.  Inputs are:

            spwidx:    spw index within file
            type:      'cross' or 'auto' (default 'cross')

        The returned array shape is (nBl/nAnt, nBin, nSpp, nPol).

        If spwidx is 'all' and the array dimensions of all spectral
        windows match, all will be returned in a single array.  In this case
        the returned dimensions will be (nBl/nAnt, nSpw, nBin, nSpp, nPol).
        """
        if type[0].lower() == 'c':
            loc = 'crossData'
        elif type[0].lower() == 'a':
            loc = 'autoData'
        else:
            raise RuntimeError('Unsupported data type')
        if spwidx == 'all':
            if not BDFSpectralWindow.dims_match(self.spws, type):
                raise RuntimeError('BDFIntegration: ' +
                                   'mixed array dimensions, spws must be ' +
                                   'retrieved indivdually')
            dshape = (-1, len(self.spws)) + self.spws[0].dshape(type)
            return self.data[loc].reshape(dshape)
        #elif ('A' in swpidx or 'C' in spwidx 
        #        or 'B' in spwidx or 'D' in spwidx):
        #    # TODO using "AC-1" syntax try to find the spw.  need
        #    # to determine whether spw.sw is swindex or sbid..
        #    # spw.swbb is name like AC_8BIT, B1D1_3BIT.
        #    pass
        spw = self.spws[spwidx]
        if loc == 'crossData':
            offs = spw.cross_offset
        elif loc == 'autoData':
            offs = spw.auto_offset
        else:
            raise RuntimeError('Unsupported data type')
        dsize = spw.dsize(type)
        dshape = (-1,) + spw.dshape(type)
        return self.data[loc][:, offs:offs+dsize].reshape(dshape)

    def get_meta(self, component, spwidx='all', corr='cross'):
        """
        Return metadata from the requested binary component (flags, actualDurations, 
        actualTimes).
        """
        # Assume the axes BAL ANT BAB SPW BIN STO
        # Note for this, npol==3 actually means 3 unlike for the data
        
        # Normalize corr input
        csize = self.numBaseline*sum([spw.msize('cross') for spw in self.spws])
        if corr[0].lower() == 'c':
            corr = 'cross'
            n0 = self.numBaseline
            offs = 0
            tsize = csize
        elif corr[0].lower() == 'a':
            corr = 'auto'
            n0 = self.numAntenna
            offs = csize
            tsize = self.numAntenna*sum([spw.msize('auto') for spw in self.spws])
        else:
            raise RuntimeError('Unsupported data type')

        data = self.data[component][offs:offs+tsize].reshape((n0,-1))

        if spwidx=='all':
            if not BDFSpectralWindow.dims_match(self.spws, corr):
                raise RuntimeError('BDFIntegration: ' +
                                   'mixed array dimensions, spws must be ' +
                                   'retrieved indivdually')
            dshape = (n0, len(self.spws)) + self.spws[0].mshape(corr)
            return data.reshape(dshape)

        spw = self.spws[spwidx]
        dshape = (n0,) + spw.mshape(corr)
        dsize = spw.msize(corr)
        if corr=='cross':
            spwoffs = spw.cross_meta_offset
        elif corr=='auto':
            spwoffs = spw.auto_meta_offset
        return data[:,spwoffs:spwoffs+dsize].reshape(dshape)

    def zerofraction(self, spwidx='all', type='cross'):
        """Returns the fraction of data points in the integration that
        are exactly zero (generally this means they have been flagged
        or otherwise not recorded by the online systems).

        Note that for WIDAR autocorrelation data, the default is to only
        record half the antennas so will typically have ~50% zeros
        according to this function."""

        if type[0].lower() == 'c':
            loc = 'crossData'
        elif type[0].lower() == 'a':
            loc = 'autoData'

        if spwidx == 'all':
            dtmp = self.data[loc].ravel()
        else:
            dtmp = self.get_data(spwidx=spwidx, type=type).ravel()
        return float(len(dtmp) - numpy.count_nonzero(dtmp))/float(len(dtmp))

# Notes for generating BDFs from scratch:
#
# Example data header:
#<sdmDataHeader xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xl="http://www.w3.org/1999/xlink" xmlns:xv="http://Alma/XVERSION" xmlns="http://Alma/XASDM/sdmbin" xsi:schemaLocation="http://Alma/XASDM/sdmbin http://almaobservatory.org/XML/XASDM/sdmbin/2/sdmDataObject.xsd" xv:schemaVersion="2" xv:revision="1.1.2.1" xv:release="ALMA-6_1_0-B" mainHeaderId="sdmDataHeader" byteOrder="Little_Endian" projectPath="0/1/1/">
#  <startTime>4922977435000000000</startTime>
#  <dataOID xl:type="locator" xl:href="uid:///evla/bdf/1416260627675" xl:title="EVLA WIDAR correlator visibility data"/>
#  <dimensionality axes="TIM">1</dimensionality>
#  <execBlock xl:href="uid:///evla/bdf/1416260627675" xl:type="simple"/>
#  <numAntenna>27</numAntenna>
#  <correlationMode>CROSS_AND_AUTO</correlationMode>
#  <spectralResolution>FULL_RESOLUTION</spectralResolution>
#  <processorType>CORRELATOR</processorType>
#  <dataStruct xsi:type="CrossAndAutoData" apc="AP_UNCORRECTED">
#    <baseband name="AC_8BIT">
#      <spectralWindow sw="1" swbb="AC_8BIT" sdPolProducts="RR RL LL" crossPolProducts="RR RL LR LL" numSpectralPoint="64" numBin="40" scaleFactor="1.000000" sideband="NOSB"/>
#    </baseband>
#    <flags size="59400" axes="BAL ANT BAB SPW BIN STO"/>
#    <actualTimes size="59400" axes="BAL ANT BAB SPW BIN STO"/>
#    <actualDurations size="59400" axes="BAL ANT BAB SPW BIN STO"/>
#    <crossData size="7188480" axes="BAL BAB SPW BIN SPP STO"/>
#    <autoData size="276480" axes="ANT BAB SPW BIN SPP STO" normalized="false"/>
#  </dataStruct>
#</sdmDataHeader>

# Build a sdmDataHeader from scratch
_nsmap_hdr = {
        'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
        'xl': 'http://www.w3.org/1999/xlink',
        'xv': 'http://Alma/XVERSION',
        None: 'http://Alma/XASDM/sdmbin'
        }


def _sdmDataHeader(time, uid, num_antenna, spws, path='0/1/1', cross=True,
                   auto=True):
    """Generate a sdmDataHeader XML element from the specified parameters:
        time: start time in SDM format (MJD ns)
        uid: unique ID for the BDF
        num_antenna: number of antennas in the data
        spws: list of BDFSpectralWindow objects; order matters!
    """
    _E = objectify.ElementMaker(annotate=False, nsmap=_nsmap_hdr)
    xl_type = '{%s}type' % _nsmap_hdr['xl']
    xsi_type = '{%s}type' % _nsmap_hdr['xsi']
    xl_href = '{%s}href' % _nsmap_hdr['xl']
    xl_title = '{%s}title' % _nsmap_hdr['xl']
    xsi_schemalocation = '{%s}schemaLocation' % _nsmap_hdr['xsi']
    xv_schemaversion = '{%s}schemaVersion' % _nsmap_hdr['xv']
    xv_revision = '{%s}revision' % _nsmap_hdr['xv']
    xv_release = '{%s}release' % _nsmap_hdr['xv']
    if cross and auto:
        corr_mode = 'CROSS_AND_AUTO'
        data_type = 'CrossAndAutoData'
    elif cross:
        corr_mode = 'CROSS_ONLY'
        data_type = 'CrossData'
    elif auto:
        corr_mode = 'AUTO_ONLY'
        data_type = 'AutoData'
    else:
        raise RuntimeError('No data type specified (cross or auto).')
    result = _E.sdmDataHeader(
            _E.startTime(time),
            _E.dataOID({
                xl_type: 'locator',
                xl_href: uid,
                xl_title: "EVLA WIDAR correlator visibility data"
                }),
            _E.dimensionality(1, axes="TIM"),
            _E.execBlock({xl_href: uid, xl_type: "simple"}),
            _E.numAntenna(num_antenna),
            _E.correlationMode(corr_mode),
            _E.spectralResolution('FULL_RESOLUTION'),
            _E.processorType('CORRELATOR'),
            _E.dataStruct({xsi_type: data_type, 'apc': 'AP_UNCORRECTED'}),
            # sdmDataHeader attributes
            {xsi_schemalocation: 'http://Alma/XASDM/sdmbin http://almaobservatory.org/XML/XASDM/sdmbin/2/sdmDataObject.xsd',
                xv_schemaversion: '2',
                xv_revision: '1.1.2.1',
                xv_release: 'ALMA-6_1_0-B',
                'mainHeaderId': 'sdmDataHeader',
                'byteOrder': 'Little_Endian',
                'projectPath': path}
            )
    # Now add the spws.  We are assuming they are sorted in the right
    # order already...
    cur_bb = None
    auto_size = 0
    cross_size = 0
    for s in spws:
        if s.swbb != cur_bb:
            cur_bb = s.swbb
            bb = _E.baseband(name=cur_bb)
            result.dataStruct.append(bb)
        bb.append(s.to_xml())
        auto_size += s.dsize('auto')
        cross_size += s.dsize('cross')
    num_baseline = num_antenna * (num_antenna-1) // 2
    auto_size *= num_antenna
    cross_size *= 2.0 * num_baseline
    if cross:
        result.dataStruct.append(_E.crossData(size='%d' % cross_size,
                                              axes="BAL BAB SPW BIN SPP STO"))
    if auto:
        result.dataStruct.append(_E.autoData(size='%d' % auto_size,
                                             axes="ANT BAB SPW BIN SPP STO",
                                             normalized="false"))
    return result

# Example data sub header:
#<sdmDataSubsetHeader xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xl="http://www.w3.org/1999/xlink" xmlns="http://Alma/XASDM/sdmbin" xsi:type="BinaryCrossAndAutoDataFXF" projectPath="0/1/1/1/">
#  <schedulePeriodTime>
#    <time>4922977454422548000</time>
#    <interval>14350776000</interval>
#  </schedulePeriodTime>
#  <dataStruct ref="sdmDataHeader"/>
#  <crossData xl:href="0/1/1/1/crossData.bin" type="FLOAT32_TYPE"/>
#  <autoData xl:href="0/1/1/1/autoData.bin"/>
#</sdmDataSubsetHeader>

# Build a sdmDataSubsetHeader element from scratch
_nsmap_subhdr = {
        'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
        'xl': 'http://www.w3.org/1999/xlink',
        None: 'http://Alma/XASDM/sdmbin'
        }


def _sdmDataSubsetHeader(time, interval, cross=True, auto=True,
                         path='0/1/1/1/'):
    _E = objectify.ElementMaker(annotate=False, nsmap=_nsmap_subhdr)
    xl_href = '{%s}href' % _nsmap_subhdr['xl']
    xsi_type = '{%s}type' % _nsmap_subhdr['xsi']
    result = _E.sdmDataSubsetHeader(
            _E.schedulePeriodTime(
                _E.time(time),
                _E.interval(interval)
                ),
            _E.dataStruct(ref='sdmDataHeader'),
            projectPath=path)
    if cross:
        result.append(_E.crossData({
            'type': 'FLOAT32_TYPE',
            xl_href: path + 'crossData.bin'
            }))
        result.attrib[xsi_type] = 'BinaryCrossDataFXF'
    if auto:
        result.append(_E.autoData({
            xl_href: path + 'autoData.bin'
            }))
        result.attrib[xsi_type] = 'BinaryAutoDataFXF'
    if cross and auto:
        result.attrib[xsi_type] = 'BinaryCrossAndAutoDataFXF'
    return result


class BDFWriter(object):
    """
    Write a BDF file.
    """
    def __init__(self, path, fname=None, bdf=None, start_mjd=None, uid=None,
                 num_antenna=None, spws=None, scan_idx=None, subscan_idx=1,
                 corr_mode=None):
        """Init BDFWrite with output filename (fname).  If the bdf
        argument contains a BDF object, its header is copied for the
        output file.  Otherwise the following arguments need to be
        filled in for the header:

            start_mjd: Start time of the BDF in MJD
            uid: UID to put in the header (eg, uid:///evla/bdf/1484080742396)
            num_antenna: Number of antennas
            spws: List of BDFSpectralWindow objects in correct order
            scan_idx: index of this scan in the SDM
            subscan_idx: idx of the subscan in the SDM (default 1)
            corr_mode:  'ca' for cross and auto, 'c' for cross, 'a' for auto
        """
        if fname is not None:
            self.fname = os.path.join(path, fname)
        else:
            self.fname = os.path.join(path,
#                                      uid.translate(string.maketrans(':/','__')))
                                      uid.replace(':/', '__').replace('/', '_'))
        self.fp = None
        self.curidx = 1
        self.mb1 = "MIME_boundary-1"
        self.mb2 = "MIME_boundary-2"
        self.len0 = 0
        self.len1 = 0
        self.len2 = 0
        self.sdmDataHeader = None
        if bdf is not None:
            self.sdmDataHeader = deepcopy(bdf.sdmDataHeader)
        else:
            cross = 'c' in corr_mode
            auto = 'a' in corr_mode
            path = '0/%d/%d/' % (scan_idx, subscan_idx)
            self.sdmDataHeader = _sdmDataHeader(int(start_mjd*86400.0e9),
                                                uid, num_antenna, spws,
                                                path=path, cross=cross,
                                                auto=auto)

    def write_header(self):
        """Open output and write the current header contents."""
        self.fp = open(self.fname, 'wb')
        tophdr = MIMEHeader()
        tophdr['MIME-Version'] = ['1.0', ]
        tophdr['Content-Type'] = ['multipart/mixed', 'boundary='+self.mb1]
        tophdr['Content-Description'] = [
                'EVLA/CORRELATOR/WIDAR/FULL_RESOLUTION', ]
        # How do we generate a new unique name?
        nsxl = self.sdmDataHeader.nsmap['xl']
        uid = self.sdmDataHeader.dataOID.attrib['{%s}href' % nsxl][5:]
        tophdr['Content-Location'] = ['http://evla.nrao.edu/wcbe/XSDM' + uid, ]
        self.fp.write(bytes(tophdr.tostring() + '\n', 'utf-8'))

        self.fp.write(bytes('--' + self.mb1 + '\n', 'utf-8'))
        xhdr = MIMEHeader()
        xhdr['Content-Type'] = ['text/xml', 'charset=utf-8']
        xhdr['Content-Location'] = ['sdmDataHeader.xml', ]
        self.fp.write(bytes(xhdr.tostring() + '\n', 'utf-8'))
        self.fp.write(etree.tostring(self.sdmDataHeader,
                                     standalone=True, encoding='utf-8') + b'\n')

    def write_integration(self, bdf_int=None, mjd=None, interval=None,
                          data=None):
        """
        Input is a BDFIntegration object (bdf_int).  The projectPath will
        be updated so that it is consistent for the file being written but
        otherwise no changes are made to the contents.

        Alternately, rather than a bdf_int, the remaining arguments can be
        filled in (for creating BDFs from scratch):

          mjd: the MJD of the midpoint of the integaration
          interval: the duration of the integration (sec)
          data: a dict whose entries are the numpy data arrays.  The
            keywords should be one or both of 'crossData' and 'autoData'
            depending on whether cross-correlations, auto-correlations, or
            both are present in the data set.
        """

        # NOTES for doing this:  bdf_int does not really need to be a
        # BDFIntegration.  It needs to act like it in the following ways:
        # 1. It needs to have a sdmDataSubsetHeader lxml Element object
        # representing the sub header.
        # 2. It needs to have a data attribute which is a dict of numpy
        # arrays containing the actual data to be written.  It

        tophdr = MIMEHeader()
        tophdr['Content-Type'] = ['multipart/related', 'boundary='+self.mb2]
        tophdr['Content-Description'] = ['data and metadata subset', ]

        ppidx = self.sdmDataHeader.attrib['projectPath'] + '%d/' % self.curidx

        hdr = MIMEHeader()
        hdr['Content-Type'] = ['text/xml', 'charset=utf-8']
        hdr['Content-Location'] = [ppidx + 'desc.xml']
        if self.len0 == 0:
            self.len0 = len(hdr.tostring()) + 12
        nxpad = self.len0 - len(hdr.tostring())
        if nxpad < 0:
            raise RuntimeError('nxpad(0)<0')
        hdr['X-pad'] = ['*'*nxpad, ]

        # Copy or generate XML sub-header
        if bdf_int is not None:
            subhdr = deepcopy(bdf_int.sdmDataSubsetHeader)
            data = bdf_int.data
        else:
            cross = 'crossData' in list(data.keys())
            auto = 'autoData' in list(data.keys())
            subhdr = _sdmDataSubsetHeader(int(mjd*86400e9), int(interval*1e9),
                                          cross=cross, auto=auto)

        subhdr.attrib['projectPath'] = ppidx
        nsxl = subhdr.nsmap['xl']
        dtypes = []
        mhdr = {}
        for dtype in ('crossData', 'autoData'):
            try:
                loc = ppidx + dtype + '.bin'
                getattr(subhdr, dtype).attrib['{%s}href' % nsxl] = loc
                dtypes += [dtype, ]
                mhdr[dtype] = MIMEHeader()
                mhdr[dtype]['Content-Type'] = ['application/octet-stream']
                mhdr[dtype]['Content-Location'] = [loc]
            except AttributeError:
                pass

        # Figure out how much X-pad to add
        subhdr_str = etree.tostring(subhdr, standalone=True, encoding='utf-8')
        if self.len1 == 0:
            self.len1 = len(subhdr_str) + len(mhdr[dtypes[0]].tostring()) + 50
        nxpad = self.len1 - (len(subhdr_str) + len(mhdr[dtypes[0]].tostring()))
        if nxpad < 0:
            raise RuntimeError('nxpad(1)<0')
        mhdr[dtypes[0]]['X-pad'] = ['*'*nxpad, ]

        # Assumes at most 2 data types.. TODO make more general?
        if len(dtypes) > 1:
            if self.len2 == 0:
                self.len2 = len(mhdr[dtypes[1]].tostring()) + 12
            nxpad = self.len2 - len(mhdr[dtypes[1]].tostring())
            if nxpad < 0:
                raise RuntimeError('nxpad(2)<0')
            mhdr[dtypes[1]]['X-pad'] = ['*'*nxpad, ]

        # TODO should check that data sizes match up with header info..

        # Now write it all out..
        self.fp.write(bytes('--' + self.mb1 + '\n', 'utf-8'))
        self.fp.write(bytes(tophdr.tostring()+'\n', 'utf-8'))

        # XML subheader part
        self.fp.write(bytes('--' + self.mb2 + '\n', 'utf-8'))
        self.fp.write(bytes(hdr.tostring() + '\n', 'utf-8'))
        self.fp.write(subhdr_str)

        # Data parts
        for dtype in dtypes:
            self.fp.write(bytes('\n--' + self.mb2 + '\n', 'utf-8'))
            self.fp.write(bytes(mhdr[dtype].tostring() + '\n', 'utf-8'))
            self.fp.write(data[dtype])

        # Close out mime
        self.fp.write(bytes('\n--' + self.mb2 + '--\n', 'utf-8'))
        self.curidx += 1

    def close(self):
        self.fp.write(bytes('--' + self.mb1 + '--\n', 'utf-8'))
        self.fp.close()
