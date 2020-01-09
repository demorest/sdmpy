from __future__ import print_function, division, absolute_import, unicode_literals # not casa compatible
from builtins import bytes, dict, object, range, map, input, int #, str # not casa compatible
from future.utils import itervalues, viewitems, iteritems, listvalues, listitems
from io import open

import os
import warnings
import numpy as np
from .scan import sdmarray

try:
    import psrchive
except ImportError:
    psrchive = None

import tempo_utils
import tempfile
def sdmpulsar_to_polyco(r,fmt='tempo_utils'):
    """Takes an SDM Pulsar row, returns a tempo_utils.polyco object."""
    p = tempo_utils.polyco(None)
    imjd = int(r.refTime/86400e9)
    fmjd = (int(r.refTime) - int(imjd*86400e9))/86400e9
    reffreq = float(r.refPulseFreq)
    refphase = float(r.refPhase)
    try:
        tspan = float(r.timeSpan)
        coeffs = sdmarray(r.phasePoly,float)
    except AttributeError:
        # Constant-period case
        tspan = 2.5*86400e9
        coeffs = [0.0,]

    p.imjd = imjd
    p.fmjd = fmjd
    p.rphase = refphase
    p.rfreq = reffreq
    p.site = '6' # assume VLA
    p.span = tspan/60e9
    p.ncoeff = len(coeffs)
    p.obsfreq = 1400.0
    p.coeffs = list(coeffs)

    if fmt=='tempo_utils':
        return p

    if fmt=='psrchive':
        if psrchive is None:
            raise RuntimeError("psrchive-format polycos requires psrchive")
        tmp = tempfile.NamedTemporaryFile(mode='w+')
        tmp.write(p.as_string())
        tmp.flush()
        pp = psrchive.polyco(tmp.name)
        tmp.close()
        return pp


# The big kludge-pile:
# _bin_file = "/lustre/aoc/sciops/pdemores/tests/phase_bins/binning_period.dat"
_bin_file = "/home/mchost/evla/scripts/opt/2016/09/16B-248/polyco/16B-248.period.log"
_bin_epochs = None   # list of MJDs
_bin_periods = None  # list of periods (s)
_mjd1970 = 40587
_clk_per_sec = int(64000000)
_clk_per_day = _clk_per_sec * int(86400)

# Routines to support pulsar binned/gated SDM data


def _get_epoch_period(mjd):
    """Get the closest binning reference epoch and period for a specified
    MJD.  This is kind of a hack until this info is stored properly in the
    SDM.  Returns a tuple (epoch, period, diff)"""
    global _bin_epochs, _bin_periods
    if _bin_epochs is None:
        _bin_epochs = []
        _bin_periods = []
        for l in open(_bin_file, 'rb'):
            # _bin_epochs += [psrchive.MJD(l.split()[0]),]
            # _bin_periods += [float(l.split()[1]),]
            if l.startswith('epoch '):
                epoch_clk = int(l.split()[1])
                mjd_int = epoch_clk//_clk_per_day + _mjd1970
                mjd_frac = (epoch_clk % _clk_per_day)/float(_clk_per_day)
            elif l.startswith('target period '):
                p = float(l.split()[3])/float(_clk_per_sec)
                _bin_epochs += [psrchive.MJD(mjd_int, mjd_frac), ]
                _bin_periods += [p, ]

    # Assumes input is a psrchive MJD.
    idx = np.argmin([abs((mjd-ep).in_seconds()) for ep in _bin_epochs])
    return (_bin_epochs[idx], _bin_periods[idx],
            (_bin_epochs[idx]-mjd).in_seconds())


# Kludge-pile version 2:
class BinLog(object):

    _mjd1970 = 40587
    _clk_per_sec = int(64000000)
    _clk_per_day = _clk_per_sec * int(86400)

    def __init__(self, sdmname,
                 logdir='/lustre/aoc/sciops/pdemores/binlogs'):
        self._logdir = logdir
        self.sdmname = sdmname
        # Following arrays contains info from each line of binlog file
        self.scanidx = []
        self.epoch = []
        self.period = []
        self._read()

    @classmethod
    def _clk_to_mjd(cls, clk):
        mjd_int = clk//cls._clk_per_day + cls._mjd1970
        mjd_frac = (clk % _clk_per_day)/float(_clk_per_day)
        return psrchive.MJD(mjd_int, mjd_frac)

    def _read(self):
        # Actually read the file
        for l in open(os.path.join(self._logdir, self.sdmname+'.binlog'), 'rb'):
            (_sdmname, _idx, _epoch_clk, _period_us, _period_clk) = l.split()
            if _sdmname != self.sdmname:
                continue
            self.scanidx.append(_idx)
            self.epoch.append(self._clk_to_mjd(int(_epoch_clk)))
            self.period.append(float(_period_clk)/float(self._clk_per_sec))

    def epoch_period(self, mjd):
        idx = np.argmin([abs((mjd-ep).in_seconds()) for ep in self.epoch])
        tdiff = (self.epoch[idx] - mjd).in_seconds()
        return self.epoch[idx], self.period[idx], tdiff


def dm_delay(dm, freq1, freq2=np.inf):
    """Return dispersion delay in seconds from freq2 to freq1 for given
    DM.  Freqs in MHz, can be inf."""
    return (dm/0.000241)*(1.0/(freq1*freq1) - 1.0/(freq2*freq2))


def rotate_phase(data, turns, axis=1, method='lin'):
    if method == 'lin':
        return rotate_phase_lin(data, turns, axis=axis)
    elif method == 'fft':
        return rotate_phase_fft(data, turns, axis=axis)


def rotate_phase_lin(data, turns, axis=1):
    """
    Rotate the data array a given number of turns, assuming the specified
    axis is pulse phase and corresponds to a full turn of phase.  Currently
    uses linear interpolation, which may be a better choice than FFT-based
    rotation for poorly resolved features.
    """

    nbin = data.shape[axis]
    r = (turns - np.floor(turns)) * float(nbin)
    ri = int(r)
    rf = r - ri
    return (1.0-rf)*np.roll(data, ri, axis=axis) + rf*np.roll(data, ri+1,
                                                              axis=axis)


def rotate_phase_fft(data, turns, axis=1):
    """
    Rotate the data array a given number of turns, assuming the specified
    axis is pulse phase and corresponds to a full turn of phase.  Uses
    FFT-based rotation.
    """

    nbin = data.shape[axis]
    shape = [1, ] * len(data.shape)
    shape[axis] = nbin
    ff = np.arange(float(nbin))
    ff[np.where(ff > nbin/2.0)] -= float(nbin)
    phs = np.exp(-2.0j*np.pi*turns*ff).reshape(shape)
    fdata = np.fft.fft(data, axis=axis)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", np.ComplexWarning)
        return np.fft.ifft(fdata*phs, axis=axis).astype(data.dtype)


def dedisperse_array(data, dm, freq, period, bin_axis=1, freq_axis=2,
                     spw_axis=None, phase_shift=0.0):
    """
    Dedisperse a generic array of data, of which one axis represents an
    entire turn of pulse phase, even sampled into bins.  freq should be an
    array giving the frequencies in MHz.  Up to two separate freq axes are
    allowed, given by the spw_axis and freq_axis arguments.  If spw_axis is
    None (not used), freq should be a 1-D array of freqs.  If spw_axis is
    set, freq should have dims (nspw, nchan).  An additional pulse phase 
    shift can be included via phase_shift (turns).
    """

    dp = -dm_delay(dm, freq)/period + phase_shift
    nchan = data.shape[freq_axis]
    fslice = [slice(None), ] * len(data.shape)
    fshape = list(data.shape)
    fshape[freq_axis] = 1
    if spw_axis is not None:
        nspw = data.shape[spw_axis]
        fshape[spw_axis] = 1
    else:
        nspw = 1
        dp = dp.reshape((1, -1))
    for ispw in range(nspw):
        if spw_axis is not None:
            fslice[spw_axis] = ispw
            dtmp0 = data.take([ispw, ], axis=spw_axis)
        for ichan in range(nchan):
            fslice[freq_axis] = ichan
            if spw_axis is None:
                dtmp = data.take(ichan, axis=freq_axis).reshape(fshape)
            else:
                dtmp = dtmp0.take([ichan, ], axis=freq_axis).reshape(fshape)
            data[tuple(fslice)] = rotate_phase(dtmp, dp[ispw, ichan],
                                        axis=bin_axis).squeeze()

# def dedisperse(scan,dm,period=None,bar=lambda x: x):
#    """Dedisperse a scan at the given DM.  If period is not specified,
#    it will be looked up in the period vs time table (not implemented)."""
#    pass
#
# def rephase(scan,parfile):
#    """Rephase a scan to a new ephemeris given in parfile."""
#    pass
