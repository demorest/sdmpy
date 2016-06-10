#! /usr/bin/env python

# calib.py -- PBD 2016/05

# Some routines to derive simple calibration solutions directly
# from data arrays.

import numpy as np
from numpy import linalg

from .bdf import ant2bl, bl2ant

def gaincal(data,axis=0,ref=0,avg=[],nit=3):
    """Derives amplitude/phase calibration factors from the data array
    for the given baseline axis.  In the returned array, the baseline
    dimension is converted to antenna.  No other axes are modified.
    Note this internally makes a transposed copy of the data so be 
    careful with memory usage in the case of large data sets.  A list
    of axes to average over before solving can be given in the avg 
    argument (length-1 dimensions are kept so that the solution can be
    applied to the original data)."""
    nbl = data.shape[axis]
    ndim = len(data.shape)
    (check,nant) = bl2ant(nbl)
    if check!=0:
        raise RuntimeError("Specified axis dimension (%d) is not a valid number of baselines" % nbl)
    for a in avg:
        data = data.mean(axis=a,keepdims=True)
    tdata = np.zeros(data.shape[:axis]+data.shape[axis+1:]+(nant,nant),
            dtype=data.dtype)
    for i in range(nbl):
        (a0,a1) = bl2ant(i)
        tdata[...,a0,a1] = data.take(i,axis=axis)
        tdata[...,a1,a0] = np.conj(data.take(i,axis=axis))
    for it in range(nit):
        (wtmp,vtmp) = linalg.eigh(tdata)
        v = vtmp[...,-1].copy()
        w = wtmp[...,-1]
        for i in range(nant):
            tdata[...,i,i] = w*(v.real[...,i]**2 + v.imag[...,i]**2)
    #result = np.sqrt(w[...,-1]).T*v[...,-1].T
    result = np.sqrt(w).T*v.T
    # First axis is now antenna.. refer all phases to reference ant
    result = (result*np.conj(result[ref])/np.abs(result[ref])).T
    # TODO try to reduce number of transposes
    outdims = range(axis) + [-1,] + range(axis,ndim-1)
    return result.transpose(outdims)

def applycal(data,caldata,axis=0,phaseonly=False):
    """Apply the complex gain calibration given in the caldata array
    to the data array.  The baseline/antenna axis must be specified in
    the axis argument.  Dimensions of all other axes must match up 
    (in the numpy broadcast sense) between the two arrays."""
    ndim = len(data.shape)
    nbl = data.shape[axis]
    nant = caldata.shape[axis]
    (check,nant_check) = bl2ant(nbl)
    if check!=0:
        raise RuntimeError("Specified axis dimension (%d) is not a valid number of baselines" % nbl)
    if nant!=nant_check:
        raise RuntimeError("Number of antennas does not match (data=%d, caldata=%d)" % (nant_check, nant))
    if phaseonly:
        caldata = caldata.copy()/abs(caldata)
        caldata[np.where(np.isfinite(caldata)==False)] = 0.0j
    # Modifies data in place.  Would it be better to return a calibrated
    # copy instead of touching the original?
    for ibl in range(nbl):
        # Must be some cleaner way to do this..?
        dslice = (slice(None),)*axis + (ibl,) + (slice(None),)*(ndim-axis-1)
        (a1,a2) = bl2ant(ibl)
        calfac = 1.0 / (caldata.take(a1,axis=axis)
                * caldata.take(a2,axis=axis).conj())
        calfac[np.where(np.isfinite(calfac)==False)] = 0.0j
        data[dslice] *= calfac


