#! /usr/bin/env python

# Average SDM/BDF data in time
# PBD 05/2016

import os, sys
import numpy as np
np.errstate(divide='ignore')
import sdmpy
import progressbar

sdmname = sys.argv[1]

sdm = sdmpy.SDM(sdmname)

tavg = 3.0 # in seconds

sdmout = sdmname + '.avg'
bdfoutpath = sdmout + '/ASDMBinary'

os.mkdir(sdmout)
os.mkdir(bdfoutpath)

for scan in sdm.scans():
    print "Processing '%s' scan %s:" % (sdmname, scan.idx)
    try:
        bdf = scan.bdf
    except IOError:
        print "Error reading bdf for scan %s, skipping" % (scan.idx,)
        continue
    # TODO maybe rename BDFs...
    # Also this assumes averaging time is the same during the whole
    # BDF.  This is always true for VLA data.
    bdfoutname = bdfoutpath + '/' + os.path.basename(scan._bdf_fname)
    navg = int(tavg / bdf[0].interval)
    nout = int(bdf.numIntegration / navg)
    delta_t = int((navg/2.0)*bdf[0].interval*1e9) # ns
    # Set up for output BDF, copying header info from the input BDF
    bdfout = sdmpy.bdf.BDFWriter(bdfoutname, bdf=bdf)
    bdfout.write_header()
    bar = progressbar.ProgressBar()
    for i in bar(range(nout)):
        bdfint = bdf[i*navg]
        count = {}
        for dtype in bdfint.data.keys():
            bdfint.data[dtype] = bdfint.data[dtype].copy()
            count[dtype] = np.zeros(bdfint.data[dtype].shape)
            count[dtype] += bdfint.data[dtype] != 0.0
        for j in range(1,navg):
            for dtype in bdfint.data.keys():
                dat = bdf[i*navg+j].data[dtype]
                bdfint.data[dtype] += dat
                count[dtype] += dat != 0.0
        for dtype in bdfint.data.keys():
            bdfint.data[dtype] /= count[dtype]
            bdfint.data[dtype][np.where(count[dtype]==0.0)] = 0.0
        # update timestamp and interval
        bdfint.sdmDataSubsetHeader.schedulePeriodTime.time += delta_t
        bdfint.sdmDataSubsetHeader.schedulePeriodTime.interval *= navg
        bdfout.write_integration(bdfint)
    bdfout.close()

    # update SDM entries with corect number of integrations, etc.
    scan._main.numIntegration = nout
    scan._main.dataSize = os.path.getsize(bdfoutname)
    scan._subscan.numIntegration = nout
    scan._subscan.numSubintegration = ('1 %d' % nout) + ' 0'*nout

# Write out new SDM tables
sdm.write(sdmout)

