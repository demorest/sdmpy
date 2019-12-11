#! /usr/bin/env python
import sys, math
import numpy as np
import logging
import sdmpy
import sdmpy.pulsar
import psrchive
#import tempo_utils

# Only needed for changing direction
try:
    import casatools
except ImportError:
    casatools = None

import argparse
par = argparse.ArgumentParser(
        description="Convert phase-binned SDM/BDF data to fold-mode PSRFITS.")
par.add_argument("sdmname", help="SDM to process")
#par.add_argument("-v", "--verbose", dest='loglevel',
#        action='store_const', const=logging.DEBUG,
#        default=logging.INFO,
#        help='Enable DEBUG output')
par.add_argument("-q", "--quiet", dest='loglevel',
        action='store_const', const=logging.WARN,
        default=logging.INFO,
        help='No INFO messages')
par.add_argument("-d", "--dm", type=float, default=0.0,
        help="dispersion measure (pc cm^-3) [%(default)s]")
par.add_argument("-p", "--period", type=float, default=0.0,
        help="period data were folded (s) [auto]")
par.add_argument("-P", "--polycos", action="store_true",
        help="use polyco file to adjust timestamps [false]")
par.add_argument("-H", "--hanning", action="store_true",
        help="apply Hanning smoothing")
par.add_argument("-D", "--direction", type=str, default='',
        help="rephase data to specified ra,dec")
args = par.parse_args()

logging.basicConfig(format="%(asctime)-15s %(levelname)8s %(message)s", 
        level=args.loglevel)

if args.direction:
    if casatools is None:
        logging.error('Rephasing sky direction requires CASA')
        sys.exit(1)
    # Could use astropy/etc to parse coords but we already 
    # require CASA..
    qa = casatools.quanta()
    ra, dec = args.direction.split(',')
    bad_coords = False
    try:
        ra_deg = qa.quantity(ra)
        dec_deg = qa.quantity(dec)
    except RuntimeError:
        bad_coords = True
    if (ra_deg['unit']!='deg' or dec_deg['unit']!='deg'
            or bad_coords):
        logging.error('Error parsing coordinates, try eg 12h34m56.7s,+12d34m56.7s')
        sys.exit(1)
    ra_rad = qa.convert(ra_deg,'rad')['value']
    dec_rad = qa.convert(dec_deg,'rad')['value']

sdmname = args.sdmname.rstrip('/')
sdm = sdmpy.SDM(sdmname, use_xsd=False)
try:
    binlog = sdmpy.pulsar.BinLog(sdmname)
except IOError:
    binlog = None

# TODO update to read 4-pol data
npol = 2

def unroll_chans(data):
    """Unroll the spw and channel dims into a single frequency axis.
    Assumes incoming dimensions are (bl, spw, bin, chan, pol)"""
    data = data.transpose((0,2,1,3,4))
    data = data.reshape((data.shape[0], 
        data.shape[1], 
        data.shape[2]*data.shape[3],
        data.shape[4]))
    return data

for scan in sdm.scans():
    if not scan.bdf.exists: continue
    if 'CALIBRATE_PHASE' in scan.intents:
        logging.info('processing cal scan %s' % scan.idx)
        # dims (bl,spw,bin,chan,pol)
        dcal = scan.bdf.get_data(scrunch=True)[...,[0,-1]].mean(2,
                keepdims=True)
        if args.hanning:
            sdmpy.calib.hanning(dcal,axis=3)
        dcal = unroll_chans(dcal)
        gcal = sdmpy.calib.gaincal(dcal,axis=0,ref=1)
    elif 'OBSERVE_TARGET' in scan.intents:
        logging.info('processing target scan %s' % scan.idx)

        psr = str(scan.source)

        # Assumes all spws have same number of bins
        nbin = int(scan.bdf.spws[0].numBin)

        # Assumes all spws have same number of chans
        freqs = scan.freqs().ravel()/1e6
        chanidx = np.argsort(freqs)
        nchan = freqs.size
        bw = sum([scan.spw(i).totBandwidth/1e6 for i in range(len(scan.spws))])

        # Initialize archive output
        arch = psrchive.Archive_new_Archive("ASP")
        arch.resize(0,npol,nchan,nbin)
        arch.set_source(psr)
        arch.set_dispersion_measure(args.dm)
        arch.set_coordinates(psrchive.sky_coord(*scan.coordinates))
        arch.set_centre_frequency(0.5*(freqs.max() + freqs.min()))
        arch.set_bandwidth(bw)
        arch.set_telescope('vla')
        # This fails in python3 for unknown reasons:
        #arch.set_state('PPQQ')
        # Workaround:
        arch.execute('e state=PPQQ')

        iout = 0

        # Try to load polycos
        try:
            if scan.pulsar is not None:
                polys = sdmpy.pulsar.sdmpulsar_to_polyco(scan.pulsar,
                        fmt='psrchive')
                # Used for testing rephasing:
                #polys0 = sdmpy.pulsar.sdmpulsar_to_polyco(scan.pulsar)
                #polys1 = tempo_utils.polycos.generate_from_polyco(
                #        "/users/pdemores/tzpar/B1937+21.par",
                #        polys0)
                logging.info('Read polycos from SDM Pulsar table')
            else:
                polycofile = '%s/%s.%d.polyco' % (binlog._logdir,
                        sdmname, int(scan.idx))
                polys = psrchive.polyco(polycofile)
                logging.info('Read polycos from %s' % polycofile)
        except Exception as ex:
            logging.info("Couldn't load polycos: " + repr(ex))
            polys = None

        bdf = scan.bdf
        arch.resize(arch.get_nsubint() + bdf.numIntegration - 1)
        for isub in range(1,bdf.numIntegration):
            logging.info("Processing subint %d/%d" % (isub,bdf.numIntegration))
            bdfsub = bdf[isub]
            dpsr = bdfsub.get_data()[...,[0,-1]]

            if args.hanning:
                sdmpy.calib.hanning(dcal,axis=3)

            dpsr = unroll_chans(dpsr)
            sdmpy.calib.applycal(dpsr,gcal,phaseonly=True)

            if args.direction:
                # Rephase data to a new sky direction
                logging.info("Rephasing sky direction")
                uvw0 = sdmpy.calib.uvw(bdfsub.time, scan.coordinates, 
                        scan.positions)
                uvw1 = sdmpy.calib.uvw(bdfsub.time, (ra_rad, dec_rad),
                        scan.positions)
                dw_us = 1e6*(uvw1[:,2] - uvw0[:,2])/299792458.0
                phs = np.outer(dw_us,freqs)
                phs = phs.reshape((phs.shape[0],1,phs.shape[1],1))
                dpsr *= np.exp(2.0j*np.pi*phs) # pos sign correct

            # Sum over baselines:
            dpsr = np.ma.masked_array(dpsr,dpsr==0.0)
            mpsr = np.real(dpsr.mean(0))

            subint = arch.get_Integration(iout)
            epoch_bdf = psrchive.MJD(float(bdfsub.time)) # only approx
            if args.period>0.0:
                #epoch = psrchive.MJD(epoch_bdf.intday())
                epoch = epoch_bdf
                p = args.period
                dt = 0.0
                logging.info('Using constant period')
            else:
                if (binlog is not None) or (polys is not None):
                    # Apply polyco time adjust
                    if args.polycos and polys is not None:
                        p = polys.period(epoch_bdf)
                        dt = (polys.period(epoch_bdf) 
                                * polys.phase(epoch_bdf).fracturns())
                        epoch = epoch_bdf - dt
                    else:
                        (epoch,p,dt) = binlog.epoch_period(epoch_bdf)
                else:
                    (epoch,p,dt) = sdmpy.pulsar._get_epoch_period(epoch_bdf)
                logging.info('Using epoch/period from dt=%.3fs' % dt)

            # These were used for testing dedispersion/rephasing:
            ## This sign is correct:
            #dphase = polys1.phase(bdfsub.time) - polys0.phase(bdfsub.time)
            #logging.warning("Testing rephasing p=%.6f dphase=%.3f"%(p,dphase))
            #sdmpy.pulsar.dedisperse_array(mpsr, args.dm, freqs, p,
            #        bin_axis=0, freq_axis=1, phase_shift=dphase)
            #sdmpy.pulsar.dedisperse_array(mpsr, args.dm, freqs_spw, p,
            #        bin_axis=2,freq_axis=3,spw_axis=1)
            #mpsr = unroll_chans(mpsr)[0,...]

            subint.set_epoch(epoch)
            subint.set_duration(bdfsub.interval)
            if p>0.0: subint.set_folding_period(p)
            if isub==0:
                wt=0.0
            else:
                wt=1.0
            for ochan in range(nchan):
                ichan = chanidx[ochan]
                subint.set_centre_frequency(ochan,freqs[ichan])
                for ipol in range(npol):
                    prof = subint.get_Profile(ipol,ochan)
                    prof.get_amps()[:] = mpsr[:,ichan,ipol]
                subint.set_weight(ochan,wt)
            iout += 1

        if polys is not None:
            arch.set_model(polys,False)
        outputname = sdmname + '.%03d.fits' % int(scan.idx)
        logging.info("unloading '%s'" % outputname)
        arch.unload(outputname)
