#! /usr/bin/env python

import os, sys
import numpy as np
import sdmpy
import sdmpy.bintab
from sdmpy.scan import sdmarray

# Convert VLA band naming to VLBA band naming
def convert_band(vla_band):
    bands = {
            'P': '90cm',
            'L': '20cm',
            'S': '13cm',
            'C': '6cm',
            'X': '4cm',
            'U': '2cm',
            'Ku': '2cm',
            'K': '1cm',
            'Q': '7mm'
            }
    try:
        vlba_band = bands[vla_band]
    except KeyError:
        vlba_band = vla_band
    return vlba_band

import argparse
par = argparse.ArgumentParser(
        description="Produce VLBA-style tsys file from a single-antenna (Y1) VLA SDM.")
par.add_argument("sdmname", help="SDM to process")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')
sdm = sdmpy.SDM(sdmname,use_xsd=False)

# Checks on number of antennas
nant = len(sdm['Antenna'])
if nant>1:
    print("Running in multi-antenna (Y3) mode", file=sys.stderr)
if (nant!=1) and (nant!=3):
    print("WARNING unexpected number of antennas (nant=%d)" % nant, 
            file=sys.stderr)

sp = sdmpy.bintab.unpacker(sdm['SysPower'])
sp.unpack()

# Note, freqs seem to go in increasing order, R, then L.
# Need to get this from the vex file generally?

#ant = 'Antenna_0' # assume one antenna
interval0 = 30.0 / 86400.0

sp_mjds = sp.row.timeMid / 86400.0e9
flags = sdm.scan(1).flags(sp_mjds,axis='ant')

if nant==1:
    ants=['Y',]
    #outf=[sys.stdout,] # old behavior
    outf=[open(sdmname+'.y.tsys','w')]
else:
    ants = []
    outf = []
    for iant in range(nant):
        ants.append('Y%d' % (iant+1))
        outf.append(open(sdmname+'.'+ants[iant].lower()+'.tsys' ,'w'))

for i in range(nant):
    print('# TSYS file created by y1swpow.py', file=outf[i])
    print('# %s %.8f %.8f' % (ants[i], sdm.scan(1).startMJD,
            sdm.scan(len(sdm['Scan'])).endMJD), file=outf[i])
    print("# antenna %s == %s" % (ants[i], str(sdm['Antenna'][i].name)), file=outf[i])
    print('# Files searched for swiched power data include:', file=outf[i])
    print('#   %s' % sdmname, file=outf[i])
    print('# ant D.O.Y. dur(days) nRecChan (tsys, bandName)[nRecChan]', file=outf[i])

for scan in sdm.scans():
    for iant in range(nant):
        print('# scan %d' % int(scan.idx), file=outf[iant])
        ant = 'Antenna_%d' % iant
        # For spws with identical freq (ie, Y3-type configs) all spws
        # get the same ID repeated 3 times here.  Taking unique values
        # avoids duplicating it in output.
        # Determine the list of unique spw IDs in increasing freq order:
        spws = [s[0] for s in sorted(set(zip(scan.spws,scan.reffreqs)),key=lambda x: x[1])]
        nspw = len(spws)
        t0 = scan.startMJD
        scandur = scan.endMJD - scan.startMJD
        navg = int(scandur/interval0)
        if navg==0: navg=1
        interval = scandur / navg
        tsys = np.zeros(nspw*2)
        bands = ['z',] * nspw * 2
        while t0 < scan.endMJD:
            for ispw, spw in enumerate(spws):
                idx = np.where((sp_mjds>t0)
                        * (sp_mjds<(t0+interval))
                        * (sp.row.antennaId==ant)
                        * (sp.row.spectralWindowId==str(spw))
                        * flags[:,iant])
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
                band = convert_band(str(sdm['Receiver'][spw].frequencyBand).replace('EVLA_',''))
                bands[ispw*2] = band
                bands[ispw*2+1] = band
            if dur>2.0/86400.: # require at least 2s of good data?
                line = '%s %.8f %.8f %d' % (ants[iant], t, dur, nspw*2)
                for ispw in range(nspw*2):
                    line += ' %.2f %s' % (tsys[ispw], bands[ispw])
                print(line, file=outf[iant])
            t0 += interval

