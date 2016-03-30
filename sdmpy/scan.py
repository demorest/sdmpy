# scan.py
#
# higher-level class to gather per-scan info from SDM and BDF sets

import string
import numpy
import os.path

from .bdf import BDF

def uid2fname(s):
    """Convert uid URL to file name (mainly for BDFs)."""
    return s.translate(string.maketrans(':/','__'))

def sdmarray(s):
    """Convert an array-valued SDM entry (string) into a numpy array."""
    fields = s.split()
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
        self._bdf_fname = os.path.join(sdm.path, 'ASDMBinary', 
            uid2fname(sdm['Main'][self.idx].dataUID))

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
    def intents(self):
        """Return the list of intents for this scan."""
        return list(sdmarray(self._scan.scanIntent))

    @property
    def antennas(self):
        """Return the list of antenna names for this scan."""
        cd = self.sdm['ConfigDescription'][self._main.configDescriptionId]
        sdm_ants = sdmarray(cd.antennaId)
        return [self.sdm['Antenna'][a].name for a in sdm_ants]

