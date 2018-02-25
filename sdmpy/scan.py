from __future__ import print_function, division, absolute_import, unicode_literals # not casa compatible
from builtins import bytes, dict, object, range, map, input#, str # not casa compatible
from future.utils import itervalues, viewitems, iteritems, listvalues, listitems
from io import open

import numpy
import os.path

from .bdf import BDF, ant2bl, bl2ant


# higher-level class to gather per-scan info from SDM and BDF sets

def uid2fname(s):
    """Convert uid URL to file name (mainly for BDFs)."""
#    return s.translate(string.maketrans(':/', '__'))
    return s.replace(':/', '__').replace('/', '_')


def sdmarray(s, dtype=None):
    """Convert an array-valued SDM entry (string) into a numpy array."""
    fields = str(s).split()
    ndim = int(fields[0])
    dims = tuple(map(int, fields[1:ndim+1]))
    return numpy.array(fields[ndim+1:], dtype=dtype).reshape(dims)


class Scan(object):
    """
    Represents a single subscan as part of a SDM/BDF dataset.  Convenience
    interface to open the BDF, get useful metadata, etc.
    """
    def __init__(self, sdm, scanidx, subscanidx=1):
        """
        sdm is the SDM object.
        scanidx is the scan number.
        subscanidx is the subscan number.
        """
        self.sdm = sdm
        self.idx = str(scanidx)
        self.subidx = str(subscanidx)
        self._bdf = None

    @property
    def bdf(self):
        if self._bdf is None:
            self._bdf = BDF(self.bdf_fname)
        return self._bdf

    @property
    def _main(self):
        """Convenience interface to the SDM Main table row."""
        return self.sdm['Main'][(self.idx, self.subidx)]

    @property
    def _scan(self):
        """Convenience interface to the SDM Scan table row."""
        return self.sdm['Scan'][self.idx]

    @property
    def _subscan(self):
        """Convenience interface to the SDM Subscan table row."""
        return self.sdm['Subscan'][(self.idx,self.subidx)]

    @property
    def _config(self):
        """Convenience interfact to the SDM ConfigDescription row."""
        return self.sdm['ConfigDescription'][self._main.configDescriptionId]

    @property
    def bdfdir(self):
        return self.sdm.bdfdir if self.sdm.bdfdir \
            else os.path.join(self.sdm.path, 'ASDMBinary')

    @property
    def bdf_fname(self):
        try:
            bdf_fname = os.path.join(self.bdfdir,
                                     uid2fname(self._main.dataUID.EntityRef.get('entityId')))
        except AttributeError:
            bdf_fname = os.path.join(self.bdfdir,
                                     uid2fname(self._main.dataOid.EntityRef.get('entityId')))
        return bdf_fname

    @property
    def source(self):
        """Source name as defined in SDM Scan table."""
        return self._scan.sourceName

    @property
    def field(self):
        """Field name as defined in SDM Subscan table."""
        return self._subscan.fieldName

    @property
    def coordinates(self):
        """
        Return the pointing coordinates (radians as given in Field table).
        """

        # Note as usual there are many redundant choices for where to
        # get this info from the SDM.  This is probably fine for standard
        # VLA observations where each scan points at a single location.
        # It may not be what is desired for OTF (mapping) type observations.
        # The SDM also seems to support some kind of polynomial in the Field
        # table; here we just return the 0th order part of this.

        return sdmarray(self.sdm['Field'][self._main.fieldId].referenceDir,
                        dtype=numpy.float)[0]

    @property
    def intents(self):
        """Return the list of intents for this scan."""
        return list(sdmarray(self._scan.scanIntent))

    @property
    def subintent(self):
        """Return the subscan intent."""
        return self._subscan.subscanIntent

    @property
    def antennas(self):
        """Return the list of antenna names for this scan."""
        sdm_ants = sdmarray(self._config.antennaId)
        return [self.sdm['Antenna'][a].name for a in sdm_ants]

    @property
    def stations(self):
        """Return the list of station names for this scan."""
        sdm_ants = sdmarray(self._config.antennaId)
        sdm_stns = [self.sdm['Antenna'][a].stationId for a in sdm_ants]
        return [self.sdm['Station'][s].name for s in sdm_stns]

    @property
    def positions(self):
        """
        Return the list of antenna posisitons (XYZ, m) for this scan.
        Result is an nant-by-3 array.
        """

        sdm_ants = sdmarray(self._config.antennaId)
        sdm_stns = [self.sdm['Antenna'][a].stationId for a in sdm_ants]
        return [sdmarray(self.sdm['Station'][s].position, dtype=numpy.float)
                for s in sdm_stns]

    @property
    def baselines(self):
        """
        Return the list of antenna pairs for this scan, in BDF ordering.
        """

        ants = self.antennas
        nant = len(ants)
        nbl = nant*(nant-1)//2
        # return ['%s-%s' % (ants[x[0]], ants[x[1]])
        #        for x in map(bl2ant, range(nbl))]
        return [(ants[x[0]], ants[x[1]])
                for x in map(bl2ant, list(range(nbl)))]

    @property
    def startMJD(self):
        return float(self._subscan.startTime/86400.0e9)

    @property
    def endMJD(self):
        return float(self._subscan.endTime/86400.0e9)

    @property
    def numIntegration(self):
        """Number of integrations as listed in the SDM Main table."""
        return int(self._main.numIntegration)

    @property
    def spws(self):
        """ Return the list of spw names """

        return [self.sdm['DataDescription'][dd_id].spectralWindowId
                for dd_id in sdmarray(self._config.dataDescriptionId)]

    @property
    def reffreqs(self):
        """ List of reference frequencies. One per spw in spws list. """

        return [float(self.spw(spwn).refFreq)
                for spwn in range(len(self.spws))]

    @property
    def numchans(self):
        """ List of number of channels per spw. One per spw in spws list. """

        return [int(self.spw(spwn).numChan) for spwn in range(len(self.spws))]

    @property
    def chanwidths(self):
        """ List of channel widths. One per spw in spws list. """

        return [float(self.spw(spwn).chanWidth)
                for spwn in range(len(self.spws))]

    def freqs(self, spwidx='all'):
        """
        Array of per-channel frequences for the given spectral window.
        If spwidx=='all', a nspw-by-nchan array will be returned giving all
        frequencies, if all spectral window have the same number of channels.
        """

        nspw = len(self.spws)
        rf = self.reffreqs
        nc = self.numchans
        cw = self.chanwidths
        if spwidx == 'all':
            if nc.count(nc[0]) != len(nc):
                raise RuntimeError("Variable number of channels")
            nc = nc[0]
            out = numpy.zeros((nspw, nc))
            for i in range(nspw):
                out[i, :] = numpy.arange(nc)*cw[i] + rf[i]
        else:
            out = numpy.arange(nc[spwidx]) * cw[spwidx] + rf[spwidx]
        return out

    def spw(self, idx):
        """Return the SpectralWindow entry for the given index in this scan."""

        dd_id = sdmarray(self._config.dataDescriptionId)[idx]
        spw_id = self.sdm['DataDescription'][dd_id].spectralWindowId
        return self.sdm['SpectralWindow'][spw_id]

    def flags(self, mjd, axis='bl', pad=False, flagval=0, expand=1.0):
        """
        Return flag array for the given time(s).  Input mjd can be scalar
        or array valued.  If axis=='ant', returned array will have dimensions
        (N_times, N_antenna), otherwise (N_times, N_baselines).  Flag array
        contains 1 for non-flagged points and flagval for flagged data.  If pad
        arg is true, four extra len-1 dimensions will be appended so that
        flags can be applied to standard bdf.get_data() results.  All
        flag start/stop times are increased by the value of the expand
        argument (seconds).
        """

        sdm_ants = self.antennas
        nant = len(sdm_ants)
        nbl = nant*(nant-1)//2
        t_ns = numpy.array(numpy.array(mjd)*86400.0e9, dtype=numpy.int64)
        exp_ns = int(expand*1e9)
        d_out = t_ns.shape
        isscalar = (d_out == ())
        if axis == 'ant':
            d_out += (nant,)
        else:
            d_out += (nbl,)
        if pad:
            d_out += (1, 1, 1, 1)
        out = numpy.ones(d_out, dtype=type(flagval))
        for flag in self.sdm['Flag']:
            tidx = numpy.where((t_ns > (int(flag.startTime)-exp_ns)) *
                               (t_ns < (int(flag.endTime)+exp_ns)))[0]
            if len(tidx) == 0:
                continue
            for a in sdmarray(flag.antennaId):
                flagant = self.sdm['Antenna'][a].name
                if flagant not in sdm_ants:
                    continue
                if axis == 'ant':
                    flagidx = [sdm_ants.index(flagant), ]
                else:
                    flagidx = numpy.where([flagant in pair
                                          for pair in self.baselines])[0]
                if isscalar:
                    out[flagidx] = flagval
                else:
                    for ii in flagidx:
                        out[tidx, ii] = flagval
        return out
