# scan.py
#
# higher-level class to gather per-scan info from SDM and BDF sets

import string
import numpy
import os.path

from .bdf import BDF, ant2bl, bl2ant

def uid2fname(s):
    """Convert uid URL to file name (mainly for BDFs)."""
    return s.translate(string.maketrans(':/','__'))

def sdmarray(s):
    """Convert an array-valued SDM entry (string) into a numpy array."""
    fields = str(s).split()
    ndim = int(fields[0])
    dims = tuple(map(int,fields[1:ndim+1]))
    return numpy.array(fields[ndim+1:]).reshape(dims)

class Scan(object):
    """
    Represents a single scan as part of a SDM/BDF dataset.  Convenience
    interface to open the BDF, get useful metadata, etc.
    """
    def __init__(self, sdm, scanidx):
        """
        sdm is the SDM object.
        scanidx is the index into the Main table.
        """
        self.sdm = sdm
        self.idx = scanidx
        self._bdf = None
        bdfdir = sdm.bdfdir if sdm.bdfdir else os.path.join(sdm.path, 'ASDMBinary')
        try:
            self._bdf_fname = os.path.join(bdfdir,
                                           uid2fname(sdm['Main'][self.idx].dataUID.EntityRef.get('entityId')))
        except AttributeError:
            self._bdf_fname = os.path.join(bdfdir,
                                           uid2fname(sdm['Main'][self.idx].dataOid.EntityRef.get('entityId')))

    @property
    def bdf(self):
        if self._bdf is None:
            self._bdf = BDF(self._bdf_fname)
        return self._bdf

    @property
    def _main(self):
        """Convenience interface to the SDM Main table row."""
        return self.sdm['Main'][self.idx]

    @property
    def _scan(self):
        """Convenience interface to the SDM Scan table row."""
        return self.sdm['Scan'][self.idx]

    @property
    def _subscan(self):
        """Convenience interface to the SDM Subscan table row."""
        return self.sdm['Subscan'][self.idx]

    @property
    def _config(self):
        """Convenience interfact to the SDM ConfigDescription row."""
        return self.sdm['ConfigDescription'][self._main.configDescriptionId]

    @property
    def source(self):
        return self._scan.sourceName

    @property
    def intents(self):
        """Return the list of intents for this scan."""
        return list(sdmarray(self._scan.scanIntent))

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
    def baselines(self):
        """Return the list of antenna pairs for this scan, in BDF ordering."""
        ants = self.antennas
        nant = len(ants)
        nbl = nant*(nant-1)/2
        #return ['%s-%s' % (ants[x[0]], ants[x[1]]) 
        #        for x in map(bl2ant, range(nbl))]
        return [(ants[x[0]], ants[x[1]]) for x in map(bl2ant, range(nbl))]

    @property
    def spws(self):
        """ Return the list of spw names """
        
        return [self.sdm['DataDescription'][dd_id].spectralWindowId for dd_id in sdmarray(self._config.dataDescriptionId)]

    @property
    def reffreqs(self):
        """ List of reference frequencies. One per spw in spws list. """

        return [self.spw(spwn).refFreq for spwn in range(len(self.spws))]

    @property
    def numchans(self):
        """ List of number of channels per spw. One per spw in spws list. """

        return [self.spw(spwn).numChan for spwn in range(len(self.spws))]

    @property
    def chanwidths(self):
        """ List of channel widths. One per spw in spws list. """

        return [self.spw(spwn).chanWidth for spwn in range(len(self.spws))]

    def spw(self,idx):
        """Return the SpectralWindow entry for the given index in this scan."""
        dd_id = sdmarray(self._config.dataDescriptionId)[idx]
        spw_id = self.sdm['DataDescription'][dd_id].spectralWindowId
        return self.sdm['SpectralWindow'][spw_id]


