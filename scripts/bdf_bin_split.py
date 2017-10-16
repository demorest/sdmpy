#! /usr/bin/env python

# Split a phase-binned SDM/BDF set into multiple single-bin
# sets.  Currently requires all scans/spws in the set to use the same
# number of bins.
# PBD 06/2016

import os, sys
from collections import namedtuple
import numpy as np
from numpy.fft import fft, ifft
np.errstate(divide='ignore')
import sdmpy
from sdmpy.pulsar import dedisperse_array
import progressbar
import argparse

par = argparse.ArgumentParser()
par.add_argument("sdmname", help="SDM to process")
par.add_argument("-d", "--dm", type=float, default=0.0,
        help="dispersion measure (pc cm^-3) [%(default)s]")
par.add_argument("-p", "--period", type=float, default=1.0,
        help="pulse period (s) [%(default)s]")
par.add_argument("-m", "--meansub", action="store_true",
        help="subtract period-averaged value from each bin")
par.add_argument("-s", "--scan", action="append", default=[],
        help="process this scan (multiple -s allowed)")
par.add_argument("-C", "--cal", action="store_true",
        help="produce sum/diff outputs (cal mode)")
par.add_argument("-A", "--fixautos", action="store_true",
        help="attempt to fix sign of autocorr poln products")
par.add_argument("-T", "--template", type=str, default='',
        help="convolve with ascii template given in filename")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')

sdm = sdmpy.SDM(sdmname)

# Get the number of bins, currently assumes it's constant throughout
# the data set.
# TODO: when processing a multi-bin data set this will need to change,
# I guess the right thing to do is get the maximum number of bins
# during the observation and set up enough output files to handle
# that case.
try:
    bdf0 = sdm.scan(args.scan[0]).bdf
except IndexError:
    bdf0 = sdm.scan(1).bdf
nbin = bdf0.spws[0].numBin

# Read a template file
tmpl = None
if args.template:
    # Assume this is a simple ascii file of template values
    tmpl = np.loadtxt(args.template,usecols=0)
    # Check that number of bins matches.  Could implement template rebinning
    # at some point.
    if len(tmpl) != nbin:
        print "Error, number of template bins (%d) does not match data (%d)" % (
                len(tmpl), nbin)
        sys.exit(1)
    # Rescale template values to preserve mean flux, and mean-subtract:
    tmpl -= tmpl.min()
    tmpl /= tmpl.sum()
    tmpl -= tmpl.mean()
    # FFT to prepare for convolution later on.  reshape so that it 
    # will broadcast correctly against a full data array.
    ftmpl = np.conj(fft(tmpl)).reshape((1,1,nbin,1,1))

# Sets up the output paths.  Use the "extra" ibin==nbin entry
# for the averaged data.
sdmout = []
bdfoutpath = []
if args.cal:
    sdmout.append(sdmname + '.sum')
    sdmout.append(sdmname + '.dif')
    nbin = 2
else:
    for ibin in range(nbin+1):
        if ibin==nbin:
            sdmout.append(sdmname + '.avg')
        else:
            sdmout.append(sdmname + '.bin%04d'%ibin)

for ibin in range(len(sdmout)):
    bdfoutpath.append(sdmout[ibin] + '/ASDMBinary')
    os.mkdir(sdmout[ibin])
    os.mkdir(bdfoutpath[ibin])

for scan in sdm.scans():
    if len(args.scan) and scan.idx not in args.scan:
        continue
    print "Processing '%s' scan %s:" % (sdmname, scan.idx)
    try:
        bdf = scan.bdf
    except IOError:
        print "Error reading bdf for scan %s, skipping" % (scan.idx,)
        continue
    if args.cal and bdf.spws[0].numBin != 2:
        print "nbin!=2 for scan %s, skipping" % (scan.idx,)
        continue
    # Array of dims (nspw,nchan) giving freqs in MHz
    freqs_spw = scan.freqs()/1e6
    bdfoutname = map(lambda x: x+'/'+os.path.basename(scan.bdf_fname), 
            bdfoutpath)
    # Set up for output BDFs, copying header info from the input BDF
    # and changing nbin to 1.
    bdfout = map(lambda x: sdmpy.bdf.BDFWriter(x,bdf=bdf), bdfoutname)
    for ibdf in bdfout: 
        for bb in ibdf.sdmDataHeader.dataStruct.baseband:
            for spw in bb.spectralWindow:
                spw.attrib['numBin'] = '1'
        # update size attributes
        for a in ('flags','actualTimes','actualDurations',
                'crossData','autoData'):
            try:
                ds = ibdf.sdmDataHeader.dataStruct.__dict__[a]
                if 'BIN' in ds.attrib['axes']:
                    sz = int(ds.attrib['size'])
                    ds.attrib['size'] = str(sz/nbin)
            except KeyError:
                pass
        ibdf.write_header()
    bar = progressbar.ProgressBar()
    for i in bar(range(bdf.numIntegration)):
        fullint = bdf[i]
        dtypes = fullint.data.keys()
        binint = []
        for ibin in range(len(bdfout)):
            # This is kind of a kludge to create a new BDFInt-like container
            # to hold the data that will be written out.  May want to clean
            # this up at some point.
            binint.append(namedtuple("NewBDFIntegration",
                ['sdmDataSubsetHeader','data'])(
                    fullint.sdmDataSubsetHeader,
                    {}
                    )
                )
            for dtype in dtypes:
                binint[-1].data[dtype] = None
        for dtype in dtypes:

            ## This version handles the general case where different
            ## spws can potentially have different numbers of channels
            ## or polarizations, etc.:
            #for ispw in range(len(fullint.spws)):
            #    # Axes should be bl/ant, bin, chan, poln
            #    data = fullint.get_data(ispw,type=dtype)
            #    for ibin in range(nbin):
            #        bindata = data.take(ibin,
            #                axis=1).reshape((data.shape[0],-1))
            #        if binint[ibin].data[dtype] is None:
            #            binint[ibin].data[dtype] = bindata.copy()
            #        else:
            #            binint[ibin].data[dtype] = np.hstack(
            #                    (binint[ibin].data[dtype], bindata) 
            #                    )

            ## This much simpler version handles the case where all spws
            ## have matching dimensions.  Will raise an error if this 
            ## is not the case.
            # Axes should be (bl/ant, spw, bin, chan, pol)
            data = fullint.get_data(type=dtype).copy()

            # extend zeros along bin axis
            dataspw = data.sum(3,keepdims=True)
            dataz = dataspw!=0.0
            dataz_m = dataspw.all(2,keepdims=True)
            data *= dataz_m
            
            # Dedisperse if needed
            if args.dm!=0.0:
                dedisperse_array(data, args.dm, freqs_spw, args.period,
                        bin_axis=2, freq_axis=3, spw_axis=1)

            # Flip sign of autocorr cross-pol imag component in 
            # every other spw.  Note I don't think there is a 
            # general way to know which spws need fixing..
            if args.fixautos:
                if ('auto' in dtype) and (data.shape[-1]==4):
                    for ispw in range(1,data.shape[1],2):
                        data[:,ispw,:,:,2] *= -1.0

            if args.cal:
                # If either one of the bins is zeroed, zero both in
                # output.
                b0 = data.take(0,axis=2)
                b1 = data.take(1,axis=2)
                wt = (b0!=0.0) * (b1!=0.0)
                binint[0].data[dtype] = (wt*(b1+b0))*0.5
                binint[1].data[dtype] = wt*(b1-b0)
            else:
                binint[nbin].data[dtype] = data.mean(axis=2)
                if tmpl is not None:
                    cdata = ifft(fft(data,axis=2)*ftmpl,axis=2)
                    for ibin in range(nbin):
                        binint[ibin].data[dtype] = \
                                cdata.take(ibin,axis=2).astype(data.dtype)
                else:
                    for ibin in range(nbin):
                        binint[ibin].data[dtype] = data.take(ibin,axis=2)
                        if args.meansub:
                            binint[ibin].data[dtype] -= binint[nbin].data[dtype]

        for ibin in range(len(bdfout)):
            bdfout[ibin].write_integration(binint[ibin])
                
    for b in bdfout: b.close()

    # update SDM entries with corect data size (should be same for all bins)
    scan._main.dataSize = os.path.getsize(bdfoutname[0])

# Write out new SDM tables
for s in sdmout: sdm.write(s)

