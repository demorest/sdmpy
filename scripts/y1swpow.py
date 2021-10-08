#! /usr/bin/env python

import os, sys
import numpy as np
import sdmpy
import sdmpy.bintab
from sdmpy.scan import sdmarray

import argparse
par = argparse.ArgumentParser(
        description="Produce VLBA-style tsys file from an single-antenna (Y1) VLA SDM.")
par.add_argument("sdmname", help="SDM to process")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')
sdm = sdmpy.SDM(sdmname,use_xsd=False)

# TODO check that it is a single antenna

sp = sdmpy.bintab.unpacker(sdm['SysPower'])
sp.unpack()

# Note, freqs seem to go in increasing order, R, then L.
# Need to get this from the vex file generally?

ant = 'Antenna_0' # assume one antenna
interval0 = 30.0 / 86400.0

sp_mjds = sp.row.timeMid / 86400.0e9
flags = sdm.scan(1).flags(sp_mjds,axis='ant')[:,0]

print('# TSYS file created by y1swpow.py')
print('# Y %.8f %.8f' % (sdm.scan(1).startMJD,
        sdm.scan(len(sdm['Scan'])).endMJD))
print('# Files searched for swiched power data include:')
print('#   %s' % sdmname)
print('# ant D.O.Y. dur(days) nRecChan (tsys, bandName)[nRecChan]')

for scan in sdm.scans():
    print('# scan %d' % int(scan.idx))
    nspw = len(scan.spws)
    t0 = scan.startMJD
    scandur = scan.endMJD - scan.startMJD
    navg = int(scandur/interval0)
    if navg==0: navg=1
    interval = scandur / navg
    tsys = np.zeros(nspw*2)
    bands = ['z',] * nspw * 2
    while t0 < scan.endMJD:
        for ispw, spw in enumerate(scan.spws):
            idx = np.where((sp_mjds>t0)
                    * (sp_mjds<(t0+interval))
                    * (sp.row.spectralWindowId==str(spw))
                    * flags)
            dur = 0.0
            if len(idx[0]):
                sp_idx = sp.row[idx]
                dur = sp_mjds[idx].max() - sp_mjds[idx].min()
                t = sp_mjds[idx].mean()
                tcal = sdmarray(sdm['CalDevice'][ant,spw].coupledNoiseCal,float)[:,0]
                tsys_idx = tcal \
                        * sp_idx.switchedPowerSum \
                        / sp_idx.switchedPowerDifference / 2.0
                tsys_idx = np.median(tsys_idx,axis=0)
                tsys[ispw*2] = tsys_idx[0]
                tsys[ispw*2+1] = tsys_idx[1]
            else:
                tsys[ispw*2] = 0.0
                tsys[ispw*2+1] = 0.0 
            band = str(sdm['Receiver'][spw].frequencyBand).replace('EVLA_','')
            bands[ispw*2] = band
            bands[ispw*2+1] = band
        if dur>2.0/86400.: # require at least 2s of good data?
            line = 'Y %.8f %.8f %d' % (t, dur, nspw*2)
            for ispw in range(nspw*2):
                line += ' %.2f %s' % (tsys[ispw], bands[ispw])
            print(line)
        t0 += interval

