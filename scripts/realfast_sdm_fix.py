#! /usr/bin/env python

# Average SDM/BDF data in time
# PBD 05/2016

import os, sys
import numpy as np
np.errstate(divide='ignore')
import lxml.objectify
import sdmpy
import progressbar

import argparse
par = argparse.ArgumentParser()
par.add_argument("sdmname", help="SDM to process")
par.add_argument("-e", "--ext", default="fix",
        help="extension to add to output SDM [%(default)s]")
par.add_argument("-b", "--bdfdir", default="",
        help="path to BDFs (optional)")
par.add_argument("-s", "--sdmonly", action="store_true",
        help="correct SDM tables only (not BDFs)")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')

sdm = sdmpy.SDM(sdmname,bdfdir=args.bdfdir,use_xsd=False)

sdmout = sdmname + '.' + args.ext
bdfoutpath = sdmout + '/ASDMBinary'

if not args.sdmonly:
    os.mkdir(sdmout)
    os.mkdir(bdfoutpath)
    for scan in sdm.scans():
        print("Processing '{0}' scan {1}:".format(sdmname, scan.idx))
        if scan.bdf.fp is None:
            print("Error reading bdf for scan {0}, skipping".format(scan.idx))
            continue

        bdf = scan.bdf

        # Set up for output BDF, copying header info from the input BDF
        bdfoutname = bdfoutpath + '/' + os.path.basename(scan.bdf_fname)
        bdfout = sdmpy.bdf.BDFWriter(bdfoutpath, fname=os.path.basename(scan.bdf_fname), bdf=bdf)
        bdfout.write_header()

        t0 = None

        bar = progressbar.ProgressBar()
        for i in bar(range(bdf.numIntegration)):
            bdfint = bdf[i]
            if i==0:
                # Get timestamp and interval info
                t0 = bdfint.sdmDataSubsetHeader.schedulePeriodTime.time
                dt = bdfint.sdmDataSubsetHeader.schedulePeriodTime.interval
            else:
                # Adjust timestamp
                bdfint.sdmDataSubsetHeader.schedulePeriodTime.time = t0 + i*dt
            # Remove all of lxml's pytype namespace cruft.  Could move this into
            # BDFWriter class by default.
            lxml.objectify.deannotate(bdfint.sdmDataSubsetHeader,
                    pytype=True, xsi=False, xsi_nil=False,
                    cleanup_namespaces=True)
            bdfout.write_integration(bdfint)
        bdfout.close()
else:
    os.mkdir(sdmout)

# Fix polarization table.  Assumes two-pol, a better way would be to
# look at what is in BDF and use that.  Also assumes circular basis.
assert len(sdm['Polarization'])==1
sdm['Polarization'][0].numCorr = 2
sdm['Polarization'][0].corrType = '1 2 RR LL'
sdm['Polarization'][0].corrProduct = '2 2 2 R R L L'

# Fix ConfigDescription to state cross-corr-only.  Not sure this
# is really needed since CASA seems OK with it either way.
assert len(sdm['ConfigDescription'])==1
sdm['ConfigDescription'][0].correlationMode = 'CROSS_ONLY'

# Write out new SDM tables
sdm.write(sdmout)

