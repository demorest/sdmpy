#! /usr/bin/env python

# Average SDM/BDF data in time
# PBD 05/2016

import os, sys
import numpy as np
np.errstate(divide='ignore')
import sdmpy
import progressbar

import argparse
par = argparse.ArgumentParser()
par.add_argument("sdmname", help="SDM to process")
par.add_argument("-t", "--time", type=float, default=3.0,
        help="averaging time (s) [%(default)s]")
par.add_argument("-e", "--ext", default="avg",
        help="extension to add to output SDM [%(default)s]")
par.add_argument("-b", "--bdfdir", default="",
        help="path to BDFs (optional)")
par.add_argument("-f", "--flagfrac", default=1.0, type=float,
        help="Flag integrations with zero fraction greater than flagfrac")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')

sdm = sdmpy.SDM(sdmname,bdfdir=args.bdfdir)

assert (args.flagfrac >= 0.) and (args.flagfrac <= 1.)

tavg = args.time # in seconds

sdmout = sdmname + '.' + args.ext
bdfoutpath = sdmout + '/ASDMBinary'

os.mkdir(sdmout)
os.mkdir(bdfoutpath)

for scan in sdm.scans():
    print("Processing '{0}' scan {1}:".format(sdmname, scan.idx))
    if scan.bdf.fp is None:
        print("Error reading bdf for scan {0}, skipping".format(scan.idx))
        continue

    bdf = scan.bdf

    # TODO maybe rename BDFs...
    # Also this assumes averaging time is the same during the whole
    # BDF.  This is always true for VLA data.
    bdfoutname = bdfoutpath + '/' + os.path.basename(scan.bdf_fname)
    navg = max(1, int(tavg / bdf[0].interval))
    nout = int(bdf.numIntegration / navg)
    delta_t = int((navg/2.0)*bdf[0].interval*1e9) # ns
    # Set up for output BDF, copying header info from the input BDF
    bdfout = sdmpy.bdf.BDFWriter(bdfoutpath, fname=os.path.basename(scan.bdf_fname), bdf=bdf)
    bdfout.write_header()
    bar = progressbar.ProgressBar()
    zeroed = 0
    for i in bar(range(nout)):
        bdfint = bdf[i*navg]
        count = {}
        for dtype in bdfint.data.keys():
            bdfint.data[dtype] = bdfint.data[dtype].copy()
            count[dtype] = np.zeros(bdfint.data[dtype].shape)
            count[dtype] += bdfint.data[dtype] != 0.0
        zerofrac = 0.
        for j in range(1, navg):
            if args.flagfrac < 1.0:
                zerofrac += bdf[i*navg+j].zerofraction()  # default looks at cross
            for dtype in bdfint.data.keys():
                dat = bdf[i*navg+j].data[dtype]
                bdfint.data[dtype] += dat
                count[dtype] += dat != 0.0
        zerofrac /= navg
        for dtype in bdfint.data.keys():
            bdfint.data[dtype] /= count[dtype]
            bdfint.data[dtype][np.where(count[dtype]==0.0)] = 0.0
            if zerofrac > args.flagfrac:
                bdfint.data[dtype][:] = 0.0  # flag data if zeros exceed limit
                zeroed += 1

        # update timestamp and interval
        bdfint.sdmDataSubsetHeader.schedulePeriodTime.time += delta_t
        bdfint.sdmDataSubsetHeader.schedulePeriodTime.interval *= navg
        bdfout.write_integration(bdfint)
    bdfout.close()

    if zeroed:
        print("Flagged {0} integrations for zerofrac exceeding flagfrac.".format(zeroed))

    # update SDM entries with corect number of integrations, etc.
    scan._main.numIntegration = nout
    scan._main.dataSize = os.path.getsize(bdfoutname)
    scan._subscan.numIntegration = nout
    scan._subscan.numSubintegration = ('1 %d' % nout) + ' 0'*nout

# Write out new SDM tables
sdm.write(sdmout)

