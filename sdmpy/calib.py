from __future__ import print_function, division, absolute_import, unicode_literals # not casa compatible
from builtins import bytes, dict, object, range, map, input#, str # not casa compatible
from future.utils import itervalues, viewitems, iteritems, listvalues, listitems
from io import open

import numpy as np
from numpy import linalg

# Use CASA for UVW calc but don't require it
try:
    import casatools
except ImportError:
    casatools = None

from .bdf import ant2bl, bl2ant

# Some routines to derive simple calibration solutions directly
# from data arrays.

def gaincal(data, axis=0, ref=0, avg=[], nit=3):
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
    (check, nant) = bl2ant(nbl)
    if check != 0:
        raise RuntimeError("Specified axis dimension (%d) is not a valid number of baselines" % nbl)
    for a in avg:
        data = data.mean(axis=a, keepdims=True)
    tdata = np.zeros(data.shape[:axis]+data.shape[axis+1:]+(nant, nant),
                     dtype=data.dtype)
    for i in range(nbl):
        (a0, a1) = bl2ant(i)
        tdata[..., a0, a1] = data.take(i, axis=axis)
        tdata[..., a1, a0] = np.conj(data.take(i, axis=axis))
    for it in range(nit):
        (wtmp, vtmp) = linalg.eigh(tdata)
        v = vtmp[..., -1].copy()
        w = wtmp[..., -1]
        for i in range(nant):
            tdata[..., i, i] = w*(v.real[..., i]**2 + v.imag[..., i]**2)
    # result = np.sqrt(w[...,-1]).T*v[...,-1].T
    result = np.sqrt(w).T*v.T
    # First axis is now antenna.. refer all phases to reference ant
    result = (result*np.conj(result[ref])/np.abs(result[ref])).T
    # TODO try to reduce number of transposes
    outdims = list(range(axis)) + [-1, ] + list(range(axis, ndim-1))
    return result.transpose(outdims)


def applycal(data, caldata, axis=0, phaseonly=False):
    """
    Apply the complex gain calibration given in the caldata array
    to the data array.  The baseline/antenna axis must be specified in
    the axis argument.  Dimensions of all other axes must match up
    (in the numpy broadcast sense) between the two arrays.
    """

    ndim = len(data.shape)
    nbl = data.shape[axis]
    nant = caldata.shape[axis]
    (check, nant_check) = bl2ant(nbl)
    if check != 0:
        raise RuntimeError("Specified axis dimension (%d) is not a valid number of baselines" % nbl)
    if nant != nant_check:
        raise RuntimeError("Number of antennas does not match (data=%d, caldata=%d)" % (nant_check, nant))
    if phaseonly:
        caldata = caldata.copy()/abs(caldata)
        caldata[np.where(np.isfinite(caldata) is False)] = 0.0j
    # Modifies data in place.  Would it be better to return a calibrated
    # copy instead of touching the original?
    for ibl in range(nbl):
        # Must be some cleaner way to do this..?
        dslice = (slice(None),)*axis + (ibl,) + (slice(None),)*(ndim-axis-1)
        (a1, a2) = bl2ant(ibl)
        calfac = 1.0 / (caldata.take(a1, axis=axis) * caldata.take(a2, axis=axis).conj())
        calfac[np.where(np.isfinite(calfac) is False)] = 0.0j
        data[dslice] *= calfac

def hanning(data, axis=0):
    """Apply hanning smoothing along the specified axis, typically this should
    be the spectral channel axis.  Modifies data array in-place."""
    data_pos = np.roll(data,1,axis=axis)
    data_neg = np.roll(data,-1,axis=axis)
    data += 0.5*(data_pos + data_neg)
    data *= 0.5

# Calculate uvw using CASA, mostly copied
# from rfpipe's calc_uvw().
def uvw(mjd, direction, antpos):
    """Return an Nbaseline-x-3 array giving U,V,W in meters for
    the given MJD, sky direction, and antenna positions.

      direction is (ra,dec) in radians
      antpos is Nant-by-3 array of antenna positions in meters.

    This can be called on an sdmpy Scan object like:

      uvw = sdmpy.calib.uvw(scan.startMJD, scan.coordinates, scan.positions)
    """
    if casatools is None:
        raise RuntimeError("")
    me = casatools.measures()
    qa = casatools.quanta()
    qq = qa.quantity
    s = me.direction('J2000',
            qq(direction[0],'rad'),
            qq(direction[1],'rad'))
    e = me.epoch('UTC', qq(mjd,'d'))
    o = me.observatory('VLA')
    me.doframe(o)
    me.doframe(e)
    me.doframe(s)
    pos = np.array(antpos)
    casapos = me.position('ITRF',
            qq(pos[:,0],'m'),
            qq(pos[:,1],'m'),
            qq(pos[:,2],'m'))
    bls = me.expand(me.touvw(me.asbaseline(casapos))[0])[1]['value']
    bls = bls.reshape((-1,3))
    # Original code from rfpipe:
    #ord1 = [(i,j) for i in range(nants) for j in range(i+1, nants)] # CASA order
    #ord2 = [(i,j) for j in range(nants) for i in range(j)]  # BDF order
    nant = pos.shape[0]
    nbl = bls.shape[0]
    # Index into the output (BDF-style) array for each value 
    # in the input (CASA-style) array.
    oidx = [ant2bl((i,j)) for i in range(nant) for j in range(i+1,nant)]
    uvw = 0.0*bls
    for i in range(nbl):
        uvw[oidx[i],:] = bls[i,:]
    return uvw


