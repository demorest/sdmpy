#! /usr/bin/env python

# Dedisperse and select individual integrations from a SDM 
# data set.
# PBD 09/2016

import os, sys, shutil
import argparse
import numpy as np
import sdmpy

par = argparse.ArgumentParser()
par.add_argument("sdmname", help="SDM to process")
par.add_argument("-i", "--idx", help="scan:int(:length) to store", 
        action="append", default=[])
par.add_argument("-e", "--ext", 
        help="extension to add for output [%(default)s]",
        default="ddcut")
par.add_argument("-d", "--dm", type=float, default=556.9,
        help="dispersion measure [%(default)s]")
par.add_argument("-t", "--time", type=float, default=0.0,
        help="keep all data with integration > time")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')

sdm = sdmpy.SDM(sdmname,use_xsd=False)

# dict of scan/integrations to keep
#keep = {'7': 9514,}
keep = {}
keeplen = {}
for i in args.idx:
    stuff = i.split(':')
    scan = stuff[0]
    iidx = stuff[1]
    try:
        length = stuff[2]
    except IndexError:
        length = 1
    keep[scan] = int(iidx)
    keeplen[scan] = int(length)

def dm_delay(dm,freq1,freq2=np.inf):
    """Return dispersion delay in seconds from freq2 to freq1 for given
    DM.  Freqs in MHz, can be inf."""
    return (dm/0.000241)*(1.0/(freq1*freq1) - 1.0/(freq2*freq2))

dm = args.dm

sdmout = sdmname + '.' + args.ext
bdfoutpath = sdmout + '/ASDMBinary'

os.mkdir(sdmout)
os.mkdir(bdfoutpath)

for scan in sdm.scans():
    print("Processing '{0}' scan {1}:".format(sdmname, scan.idx))
    try:
        bdf = scan.bdf
    except IOError:
        print("Error reading bdf for scan {0}, skipping".format(scan.idx))
        continue

    bdfoutname = bdfoutpath + '/' + os.path.basename(scan.bdf_fname)

    # Check either intent or int time, copy full scan if matching
    if (args.time>0 and bdf[0].interval>args.time):
        print("  copying scan {0}".format(scan.idx))
        shutil.copy(scan.bdf_fname, bdfoutname)
        continue

    if scan.idx not in keep.keys():
        print("  no integrations to keep for scan {0}, skipping"
              .format(scan.idx))
        # Fill in X1
        scan._main.dataUID.EntityRef.attrib['entityId'] = b'uid://evla/bdf/X1'
        continue

    # If we got here, need to select appropriate integration(s), dedisperse
    # and output.

    dt = bdf[0].interval # in sec
    freqs = scan.freqs()/1e6
    delays = dm_delay(dm,freqs,freqs.max())
    delays_samp = np.rint(delays/dt)
    uniq_delays_samp = sorted(set(delays_samp.ravel()))
    nint_read = delays_samp.max()+1 # number of integrations to read per event

    nout = keeplen[scan.idx]
    for offs in range(nout):

        int0 = keep[scan.idx] + offs
        bdfint = bdf[int0]

        int_time = bdfint.sdmDataSubsetHeader.schedulePeriodTime.time
        int_interval = bdfint.sdmDataSubsetHeader.schedulePeriodTime.interval

        if offs==0:
            t0 = int_time - int(int_interval)//2
            bdfout = sdmpy.bdf.BDFWriter(bdfoutpath,
                                         fname=os.path.basename(scan.bdf_fname),
                                         bdf=bdf)
            bdfout.sdmDataHeader.startTime = t0
            bdfout.write_header()

        t1 = int_time + int(int_interval)//2

        for dtype in bdfint.data.keys():
            # Copies the zero-delay data array
            bdfint.data[dtype] = bdfint.data[dtype].copy()
            dat0 = bdfint.get_data(type=dtype)
            # loop over the delays, copy data
            for ii in uniq_delays_samp:
                if (int0+ii)>=bdf.numIntegration: continue
                dat = bdf[int(int0+ii)].get_data(type=dtype)
                fidx = np.where(delays_samp == ii)
                # May be a more clever numpy-ish way to do this?
                for jj in range(len(fidx[0])):
                    dat0[:,fidx[0][jj],0,fidx[1][jj],:] = \
                            dat[:,fidx[0][jj],0,fidx[1][jj],:]

        # Write this integration
        bdfout.write_integration(bdfint)

    bdfout.close()

    # update SDM entries with corect number of integrations
    scan._main.numIntegration = nout
    scan._main.dataSize = os.path.getsize(bdfoutname)
    scan._subscan.numIntegration = nout
    scan._subscan.numSubintegration = ('1 %d' % nout) + ' 0'*nout

    # Why must time be stored in four separate places...??
    scan._main.time = t0  # correct?
    scan._main.interval = int_interval*nout
    scan._scan.startTime = t0
    scan._scan.endTime = t1
    scan._subscan.startTime = t0
    scan._subscan.endTime = t1

# Write out new SDM tables
sdm.write(sdmout)
