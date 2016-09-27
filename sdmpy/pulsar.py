#! /usr/bin/env python

# pulsar.py -- PBD 2016/08

# Routines to support pulsar binned/gated SDM data

import numpy as np

# TODO maybe reconsider this dependency, although the MJD class is convenient
from psrchive import MJD

# The big kludge-pile:
_bin_file = "/lustre/aoc/sciops/pdemores/tests/phase_bins/binning_period.dat"
_bin_epochs = None   # list of MJDs
_bin_periods = None  # list of periods (s)

def _get_epoch_period(mjd):
    """Get the closest binning reference epoch and period for a specified
    MJD.  This is kind of a hack until this info is stored properly in the 
    SDM.  Returns a tuple (epoch, period, diff)"""
    if _bin_epochs is None:
        for l in open(_bin_file):
            _bin_epochs += [MJD(l.split()[0]),]
            _bin_periods += [float(l.split()[1]),]
    # Assumes input is a psrchive MJD.  
    idx = np.argmin([abs((mjd-ep).in_seconds()) for ep in _bin_epochs])
    return (_bin_epochs[idx], _bin_periods[idx], 
            (_bin_epochs[idx]-mjd).in_seconds())

def dm_delay(dm,freq1,freq2=np.inf):
    """Return dispersion delay in seconds from freq2 to freq1 for given
    DM.  Freqs in MHz, can be inf."""
    return (dm/0.000241)*(1.0/(freq1*freq1) - 1.0/(freq2*freq2))

def rotate_phase(data,turns,axis=1):
    """Rotate the data array a given number of turns, assuming the specified
    axis is pulse phase and corresponds to a full turn of phase.  Currently
    uses linear interpolation, which may be a better choice than FFT-based
    rotation for poorly resolved features."""
    nbin = data.shape[axis]
    r = (turns - np.floor(turns)) * float(nbin)
    ri = int(r)
    rf = r - ri
    return (1.0-rf)*np.roll(data,ri,axis=axis)+rf*np.roll(data,ri+1,axis=axis)

def dedisperse_array(data,dm,freq,period,bin_axis=1,freq_axis=2,spw_axis=None):
    """Dedisperse a generic array of data, of which one axis represents an
    entire turn of pulse phase, even sampled into bins.  freq should be an
    array giving the frequencies in MHz.  Up to two separate freq axes are
    allowed, given by the spw_axis and freq_axis arguments.  If spw_axis is 
    None (not used), freq should be a 1-D array of freqs.  If spw_axis is 
    set, freq should have dims (nspw, nchan)."""
    dp = -dm_delay(dm,freq)/period
    nchan = data.shape[freq_axis]
    fslice = [slice(None),] * len(data.shape)
    fshape = list(data.shape)
    fshape[freq_axis] = 1
    if spw_axis is not None:
        nspw = data.shape[spw_axis]
        fshape[spw_axis] = 1
    else:
        nspw = 1
        dp = dp.reshape((1,-1))
    for ispw in range(nspw):
        if spw_axis is not None: 
            fslice[spw_axis] = ispw
            dtmp0 = data.take([ispw,],axis=spw_axis)
        for ichan in range(nchan):
            fslice[freq_axis] = ichan
            if spw_axis is None:
                dtmp = data.take(ichan,axis=freq_axis).reshape(fshape)
            else:
                dtmp = dtmp0.take([ichan,],axis=freq_axis).reshape(fshape)
            data[fslice] = rotate_phase(dtmp,dp[ispw,ichan],
                    axis=bin_axis).squeeze()

#def dedisperse(scan,dm,period=None,bar=lambda x: x):
#    """Dedisperse a scan at the given DM.  If period is not specified,
#    it will be looked up in the period vs time table (not implemented)."""
#    pass
#
#def rephase(scan,parfile):
#    """Rephase a scan to a new ephemeris given in parfile."""
#    pass
