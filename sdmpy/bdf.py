#! /usr/bin/env python

# bdf.py -- PBD 2015/06

# This class provides access to EVLA/ALMA Binary Data Format (BDF)
# files.  The approach used here is based loosely on 'bdfparse.py' 
# originally by Peter Williams.  So far this is somewhat targeted
# to EVLA/WIDAR data files, it may not handle all the ALMA variants.

try:
    from lxml import etree
except ImportError:
    from xml.etree import ElementTree as etree

import os
import sys
import string
import re
import mmap
import numpy
from collections import namedtuple

class MIMEPart(namedtuple('MIMEPart','hdr body')):
    """
    Simple class for representing one part of a MIME message.
    Has two member variable:

      hdr  = Dict of MIME header key/value pairs
      body = Body of message.  In our usage, can be a file offset in bytes
             (for binary parts), a string (for text) or a list of MIMEPart
             objects (for multipart).

    The loc property is a shortcut for the Content-Location header
    parameter.

    The type property is a shortcut for Content-Type
    """

    @property
    def loc(self):
        try:
            return self.hdr['Content-Location'][0]
        except KeyError:
            return None

    @property
    def type(self):
        try:
            return self.hdr['Content-Type'][0]
        except KeyError:
            return None

def basename_noext(path):
    return os.path.basename(os.path.splitext(path)[0])

# TODO find a better way to get the namespace automatically?
_ns = '{http://Alma/XASDM/sdmbin}'
def _stripns(tag):
    return re.sub('{.+}','',tag)


class BDF(object):
    """
    Class representing a single BDF file.  For example:

        b = bdf.BDF('uid____evla_bdf_1433189755525')

    Individual integration data is returned as BDFIntegration objects via
    either b.get_integration(idx) or b[idx].  Other useful methods include:

        b.basebands      # list of baseband ids
        b.spws           # dict of spectral windows per baseband
        b.numAntenna     # number of antennas
        b.numBaseline    # number of baselines
        b.numIntegration # number of integrations in file
        b.sdmDataHeader  # XML etree version of full header

    """

    def __init__(self, fname):
        self.fname = fname
        self.fp = open(fname, 'r')
        self.mmdata = mmap.mmap(self.fp.fileno(), 0, mmap.MAP_PRIVATE,
                mmap.PROT_READ)
        self.read_mime()
        self.parse_spws()

    @staticmethod
    def split_mime(line):
        idx = line.index(':')
        key = line[:idx]
        vals = map(string.strip, line[idx+1:].split(';'))
        return (key, vals)

    @staticmethod
    def mime_boundary(mime_hdr):
        if mime_hdr['Content-Type'][0].startswith('multipart/'):
            for v in mime_hdr['Content-Type']:
                if v.startswith('boundary='):
                    return v[v.index('=')+1:]
        return None

    def read_mime_part(self,boundary=None,recurse=False):
        """
        Read a MIME content part starting at the current file location.
        Return value is a MIMEPart object, which has elements:

            hdr    dict of MIME header key / value pairs

            body   string if Content-Type was 'text/xml', offset into
                   the file if 'application/octet-stream', or list of
                   other MIMEParts for a 'multipart/*'.

        If recurse is True, will read/return the contents of a multipart
        (and any multiparts found at lower levels).  Otherwise will read
        one header/body unit and pointer will be left at the start of 
        the next one (or first sub-part for multiparts).
        """
        hdr = {}
        body = None
        in_hdr = True
        binary_type = False
        multipart_type = False
        # Note, need to use readline() rather than iterating over file
        # because we need to recover file positions and seek ahead.
        # The "for line in file" type loop reads ahead so is not compatible
        # with this approach.
        while True:

            line = self.fp.readline()

            # hit EOF
            if line=='':
                return MIMEPart(hdr, body)

            # Check for multipart boundary marker
            if boundary is not None:
                if in_hdr:
                    # If we are starting, ignore a 'start' marker,
                    # quit on a 'done' marker
                    if line=='--'+boundary+'\n':
                        continue
                    elif line=='--'+boundary+'--\n':
                        return MIMEPart({}, None)
                else:
                    # This marks the end of a part, rewind so that the 
                    # next part can be parsed, and return results
                    if line.startswith('--' + boundary):
                        self.fp.seek(-len(line),1)
                        return MIMEPart(hdr, body)

            if line=='\n':
                # Got blank line, the next part will be body.  We
                # want to skip it if this is a binary part, otherwise
                # read and return the body.
                in_hdr = False
                if binary_type:
                    # Note the location within the file and skip
                    # ahead by the correct amount.
                    bin_name = basename_noext(hdr['Content-Location'][0])
                    body = self.fp.tell()
                    # Need to add one extra byte for the newline
                    self.fp.seek(self.bin_size[bin_name]+1, 1)
                elif multipart_type:
                    if recurse:
                        # Parse the parts and add to a list
                        while True:
                            #print "recur b='%s'" % boundary 
                            pmime = self.read_mime_part(boundary=boundary,
                                        recurse=True)
                            #print phdr
                            #print pbody
                            if pmime.hdr == {}:
                                return MIMEPart(hdr,body)
                            else:
                                body.append(pmime)
                continue

            if in_hdr:
                # Still in the header, parse the line as MIME key/val
                (key, vals) = self.split_mime(line)
                hdr[key] = vals
                if key=='Content-Type':
                    if vals[0].startswith('multipart/'):
                        multipart_type = True
                        boundary = self.mime_boundary(hdr)
                        body = []
                    elif vals[0] == 'application/octet-stream':
                        binary_type = True
            else:
                if not binary_type:
                    # In body part of a non-binary type
                    if body is None: body = line
                    else: body += line
                else:
                    # Should not really get here, means size calculation
                    # failed or file is otherwise messed up... what to do?
                    raise RuntimeError('BDF MIME parsing failure')

    # Size in bytes of each data element. In principle, the crossData 
    # type needs to be read from the headers while all others have 
    # pre-set values...
    bin_dtype_size = {
            'flags':           4, # INT32
            'actualTimes':     8, # INT64
            'actualDurations': 8, # INT64
            'zeroLags':        4, # FLOAT32
            'autoData':        4, # FLOAT32
            'crossData':       4, # FLOAT32 but should be determined from file
            }

    # Basic data type for each array.  Note, does not necessarily
    # match up with the sizes above due to how BDF defines the 'size'
    # attribute..
    bin_dtype = {
            'autoData':        numpy.float32,
            'crossData':       numpy.complex64,
            }

    def read_mime(self):
        self.fp.seek(0,0) # Go back to start
        if not self.fp.readline().startswith('MIME-Version:'):
            raise RuntimeError('Invalid BDF: missing MIME-Version')

        # First we need to read and parse only the main XML header in order
        # to get sizes of the binary parts.  Note, the info stored in 
        # self.bin_size is in bytes, rather than the weird BDF units.
        mime_hdr = self.read_mime_part().hdr
        sdmDataMime = self.read_mime_part(boundary=self.mime_boundary(mime_hdr))
        if sdmDataMime.loc != 'sdmDataHeader.xml':
            raise RuntimeError('Invalid BDF: missing sdmDataHeader.xml')
        self.sdmDataHeader = etree.fromstring(sdmDataMime.body)
        self.bin_size = {}
        self.bin_axes = {}
        for e in self.sdmDataHeader.iter():
            if 'size' in e.attrib.keys() and 'axes' in e.attrib.keys():
                binname = _stripns(e.tag)
                self.bin_size[binname] = int(e.attrib['size']) \
                        * self.bin_dtype_size[binname]
                self.bin_axes[binname] = e.attrib['axes'].split()

        # Then we go back to the beginning and parse the whole MIME
        # structure.
        self.fp.seek(0,0) # reset to start again
        full_mime = self.read_mime_part(recurse=True)
        self.mime_ints = full_mime.body[1:]

    @property
    def projectPath(self):
        return self.sdmDataHeader.attrib['projectPath']

    @property
    def numIntegration(self):
        return len(self.mime_ints)

    @property
    def numAntenna(self):
        return int(self.sdmDataHeader.find(_ns+'numAntenna').text)

    @property
    def numBaseline(self):
        return (self.numAntenna*(self.numAntenna-1))/2

    def parse_spws(self):
        self.basebands = []
        self.spws = {}
        dataStruct = self.sdmDataHeader.find(_ns+'dataStruct')
        # Offsets into the cross and data arrays where this spw
        # can be found.  This info could be reconstructed later
        # but seems convenient to do it here.
        cross_offset = 0
        auto_offset = 0
        for bb in dataStruct.iterfind(_ns+'baseband'):
            bbname = bb.attrib['name']
            self.basebands.append(bbname)
            self.spws[bbname] = []
            for spw_elem in bb.iterfind(_ns+'spectralWindow'):
                # Build a list of spectral windows for each baseband
                spw = SpectralWindow(spw_elem,cross_offset,auto_offset)
                cross_offset += spw.dsize('cross')
                auto_offset += spw.dsize('auto')
                self.spws[bbname].append(spw)

    def get_integration(self,idx):
        return BDFIntegration(self,idx)

    def __getitem__(self,idx):
        return self.get_integration(idx)

class SpectralWindow(object):
    """Spectral window class.  Initialize from the XML element."""

    def __init__(self, spw_elem, cross_offset=None, auto_offset=None):
        self._attrib = spw_elem.attrib
        self.cross_offset = cross_offset
        self.auto_offset = auto_offset

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

    def npol(self,type):
        """Return number of polarization array elements for the given data
        type (cross or auto)."""
        try:
            if type[0].lower()=='c':
                return len(self._attrib['crossPolProducts'].split())
            elif type[0].lower()=='a':
                # Good enough?
                l = len(self._attrib['sdPolProducts'].split())
                return 4 if l==3 else l # 3==4 in BDF math! :)
        except KeyError:
            return 0

    def dshape(self,type):
        """Return shape tuple of data array for this spectral window, 
        in number of data elements (real for auto, complex for cross)."""
        return (self.numBin, self.numSpectralPoint, self.npol(type))

    def dsize(self,type):
        """Return size of data array for this spectral window, in number of
        data elements (real for auto, complex for cross)."""
        return numpy.product(self.dshape(type))

class BDFIntegration(object):
    """
    Describes and holds data for a single intgration within a BDF file.
    This should be derived from an existing BDF object using 
    get_integration() or indexing, ie:

        b = bdf.BDF('some_file')

        # Get the 5th integration, these two are equivalent:
        i = b.get_integration(5)
        i = b[5]

        # Read the cross-corr data array for spectral window 0 in 
        # the AC baseband:
        dat = i.get_data('AC_8BIT',0)

    Other potentially useful info:

        i.basebands            # list of baseband IDs
        i.spws                 # dict of spws per baseband
        i.numAntenna           # obvious
        i.numBaseline          # "
        i.sdmDataSubsetHeader  # XML etree version of full sub-header

    """

    def __init__(self,bdf,idx):
        # Get the main header
        self.sdmDataSubsetHeader = etree.fromstring(
                bdf.mime_ints[idx].body[0].body)
        # Copy some info from the BDF headers
        self.basebands = bdf.basebands
        self.spws = bdf.spws
        self.bin_axes = bdf.bin_axes
        self.numAntenna = bdf.numAntenna
        self.numBaseline = bdf.numBaseline
        # Get the binary data
        self.data = {}
        for m in bdf.mime_ints[idx].body[1:]:
            btype = basename_noext(m.loc)
            bsize = bdf.bin_size[btype] # size of the binary blob in bytes
            baxes = self.bin_axes[btype]
            # Determine outer step size for the array, either baselines
            # antennas of baseline+antenna.  We can't apply the other
            # dimensions here because the number of elements can vary
            # per spw.
            if baxes[0]=='BAL' and baxes[1]=='ANT':
                shape=(self.numBaseline+self.numAntenna,-1)
            elif baxes[0]=='BAL':
                shape=(self.numBaseline,-1)
            elif baxes[0]=='ANT':
                shape=(self.numAntenna,-1)
            else:
                shape=(-1,) # Don't know what to do, just leave flat array
            self.data[m.loc] = numpy.fromstring(bdf.mmdata[m.body:m.body+bsize],
                    dtype=bdf.bin_dtype[btype]).reshape(shape)

    @property
    def projectPath(self):
        return self.sdmDataSubsetHeader.attrib['projectPath']

    def get_data(self,baseband,spwidx,type='cross'):
        """
        Return the data array for the given subset.  Inputs are:

            baseband:  baseband ID string
            spwidx:    spw index within baseband
            type:      'cross' or 'auto' (default 'cross')

        The returned array shape is (nBl/nAnt, nBin, nSpp, nPol).
        """
        spw = self.spws[baseband][spwidx]
        if type[0].lower()=='c': 
            loc = self.projectPath + 'crossData.bin'
            offs = spw.cross_offset
        elif type[0].lower()=='a': 
            loc = self.projectPath + 'autoData.bin'
            offs = spw.auto_offset
        else:
            raise RuntimeError('Unsupported data type')
        dsize = spw.dsize(type)
        dshape = (-1,) + spw.dshape(type)
        return self.data[loc][:,offs:offs+dsize].reshape(dshape)

