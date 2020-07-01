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
import tempo_utils
from sdmpy.pulsar import dedisperse_array, sdmpulsar_to_polyco
import sdmpy.calib
import progressbar
import argparse

par = argparse.ArgumentParser()
par.add_argument("sdmname", help="SDM to process")
par.add_argument("-d", "--dm", type=float, default=0.0,
        help="dispersion measure (pc cm^-3) [%(default)s]")
par.add_argument("-p", "--period", type=float, default=None,
        help="override pulse period (s) [%(default)s]")
par.add_argument("-m", "--meansub", action="store_true",
        help="subtract period-averaged value from each bin")
par.add_argument("-e", "--extend", action="store_true",
        help="extend 1-bin scans across all output bins")
par.add_argument("-s", "--scan", action="append", default=[],
        help="process this scan (multiple -s allowed)")
par.add_argument("-C", "--cal", action="store_true",
        help="produce sum/diff outputs (cal mode)")
par.add_argument("-A", "--fixautos", action="store_true",
        help="attempt to fix sign of autocorr poln products")
par.add_argument("-T", "--template", type=str, default='',
        help="convolve with ascii template given in filename")
par.add_argument("-E", "--ephemeris", type=str, default='',
        help="rephase using specified parfile")
par.add_argument("-H", "--hanning", action="store_true",
        help="Hanning-smooth before dedispersion")
par.add_argument("-M", "--memory", type=float, default=sdmpy.bdf._mmap_limit/(1<<30),
        help="Memory buffer size for BDFs (GB) [%(default)s]")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')

sdm = sdmpy.SDM(sdmname,use_xsd=False)

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
#nbin = bdf0.spws[0].numBin
nbin = max([s.bdf.spws[0].numBin for s in sdm.scans() if s.bdf.exists])

sdmpy.bdf._mmap_limit = int(args.memory * (1<<30))

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

# Sets up the output paths.  bin 0 is for the averaged data.
sdmout = []
bdfoutpath = []
if args.cal:
    sdmout.append(sdmname + '.sum')
    sdmout.append(sdmname + '.dif')
    nbin = 2
else:
    for ibin in range(nbin+1):
        if ibin==0:
            sdmout.append(sdmname + '.avg')
        else:
            sdmout.append(sdmname + '.bin%04d'%(ibin-1))

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
    if not scan.bdf.exists:
        print "Missing bdf for scan %s, skipping" % (scan.idx,)
        continue
    if args.cal and bdf.spws[0].numBin != 2:
        print "nbin!=2 for scan %s, skipping" % (scan.idx,)
        continue

    # Get number of bins in scan, set up output nbin
    nbin_scan = bdf.spws[0].numBin
    if args.extend and nbin_scan==1:
        nbin_out = nbin
    else:
        nbin_out = nbin_scan

    # Array of dims (nspw,nchan) giving freqs in MHz
    freqs_spw = scan.freqs()/1e6

    # Get the polycos as needed
    do_rephase = False
    polys = None
    if scan.pulsar is not None:
        polys = sdmpulsar_to_polyco(scan.pulsar)
        if args.ephemeris:
            polys_new = tempo_utils.polycos.generate_from_polyco(
                    args.ephemeris, polys)
            if len(polys_new)==0:
                print "Error generating polycos from '%s'" % args.ephemeris
                print "Tempo output:"
                print polys_new.tempo_output
                sys.exit(1)
            do_rephase = True
    if (nbin_scan>1 and args.dm!=0.0 and 
            polys is None and args.period is None):
        print "Missing pulse period info for scan %s, skipping" % (scan.idx,)
        continue

    # Set up for output BDFs, copying header info from the input BDF
    # and changing nbin to 1.
    bdfoutname = map(lambda x: x+'/'+os.path.basename(scan.bdf_fname), 
            bdfoutpath[:nbin_out+1])
    bdfout = map(lambda x: sdmpy.bdf.BDFWriter('',x,bdf=bdf), bdfoutname)
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
                    ds.attrib['size'] = str(sz/nbin_scan)
            except KeyError:
                pass
        ibdf.write_header()

    bar = progressbar.ProgressBar()
    for i in bar(range(bdf.numIntegration)):
        fullint = bdf[i]
        dtypes = fullint.data.keys()
        binint = []
        if args.period is not None:
            period = args.period
        elif polys is not None:
            period = 1.0/polys.freq(fullint.time)
        else:
            period = 0.0
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

            # Hanning-smooth if requested
            if args.hanning:
                sdmpy.calib.hanning(data,axis=3)

            # extend zeros along bin axis
            dataspw = data.sum(3,keepdims=True)
            dataz = dataspw!=0.0
            dataz_m = dataspw.all(2,keepdims=True)
            data *= dataz_m

            # get phase shift if rephasing
            if do_rephase:
                dphase = polys_new.phase(fullint.time) \
                        - polys.phase(fullint.time)
            else:
                dphase = 0.0
            
            # Dedisperse if needed
            if args.dm!=0.0 and nbin_scan>1:
                dedisperse_array(data, args.dm, freqs_spw, period,
                        bin_axis=2, freq_axis=3, spw_axis=1,
                        phase_shift=dphase)

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
                binint[0].data[dtype] = data.mean(axis=2)
                if tmpl is not None:
                    cdata = ifft(fft(data,axis=2)*ftmpl,axis=2)
                    for ibin in range(nbin_scan):
                        binint[ibin+1].data[dtype] = \
                                cdata.take(ibin,axis=2).astype(data.dtype)
                else:
                    for ibin in range(nbin_scan):
                        binint[ibin+1].data[dtype] = data.take(ibin,axis=2)
                        if args.meansub:
                            binint[ibin+1].data[dtype] -= binint[0].data[dtype]

        for ibin in range(len(bdfout)):
            if args.extend and nbin_scan==1:
                bdfout[ibin].write_integration(binint[0])
            else:
                bdfout[ibin].write_integration(binint[ibin])
                
    for b in bdfout: b.close()

    # update SDM entries with corect data size (should be same for all bins)
    scan._main.dataSize = os.path.getsize(bdfoutname[0])

# Write out new SDM tables
for s in sdmout: sdm.write(s)

