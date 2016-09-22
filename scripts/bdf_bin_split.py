#! /usr/bin/env python

# Split a phase-binned SDM/BDF set into multiple single-bin
# sets.  Currently requires all scans/spws in the set to use the same
# number of bins.
# PBD 06/2016

import os, sys
from collections import namedtuple
import numpy as np
np.errstate(divide='ignore')
import sdmpy
import progressbar

sdmname = sys.argv[1]

sdm = sdmpy.SDM(sdmname)

# Get the number of bins, currently assumes it's constant throughout
# the data set.
bdf0 = sdm.scan(1).bdf
nbin = bdf0.spws[0].numBin

# Sets up the output paths.  Use the "extra" ibin==nbin entry
# for the averaged data.
sdmout = []
bdfoutpath = []
for ibin in range(nbin+1):
    if ibin==nbin:
        sdmout.append(sdmname + '.avg')
    else:
        sdmout.append(sdmname + '.bin%04d'%ibin)
    bdfoutpath.append(sdmout[ibin] + '/ASDMBinary')
    os.mkdir(sdmout[ibin])
    os.mkdir(bdfoutpath[ibin])

for scan in sdm.scans():
    print "Processing '%s' scan %s:" % (sdmname, scan.idx)
    try:
        bdf = scan.bdf
    except IOError:
        print "Error reading bdf for scan %s, skipping" % (scan.idx,)
        continue
    bdfoutname = map(lambda x: x+'/'+os.path.basename(scan._bdf_fname), 
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
        for ibin in range(nbin+1):
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
            data = fullint.get_data(type=dtype)
            for ibin in range(nbin):
                binint[ibin].data[dtype] = data.take(ibin,axis=2)
            binint[nbin].data[dtype] = data.mean(axis=2)

        for ibin in range(nbin+1):
            bdfout[ibin].write_integration(binint[ibin])
                
    for b in bdfout: b.close()

    # update SDM entries with corect data size (should be same for all bins)
    scan._main.dataSize = os.path.getsize(bdfoutname[0])

# Write out new SDM tables
for s in sdmout: sdm.write(s)

