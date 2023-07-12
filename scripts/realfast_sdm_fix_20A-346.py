#! /usr/bin/env python
import os, sys
import sdmpy

import argparse
par = argparse.ArgumentParser()
par.add_argument("sdmname", help="SDM to process")
par.add_argument("-e", "--ext", default="fix",
        help="extension to add to output SDM [%(default)s]")
args = par.parse_args()

sdmname = args.sdmname.rstrip('/')
sdmout = sdmname + '.' + args.ext

print("Processing '%s'" % sdmname)

if '20A-346' not in sdmname:
    # Could try harder to protect users from themselves but probably not worth it
    print("WARNING: This script is intended for realfast SDMs from project 20A-346 only")

sdm = sdmpy.SDM(sdmname, use_xsd=False)

# List of spws to keep using SDM table numbering.  Could try to determine
# this from BDF header, don't know how reliable that is.
rf_spws = [21, 22, 23, 27, 28, 29, 30, 31]
spwid_keep = ['SpectralWindow_%d'%i for i in rf_spws]

def clean_rows(table,keep,idtag='spectralWindowId'):
    to_remove = [r for r in table if str(r.__getattr__(idtag)) not in keep]
    for r in to_remove: r.getparent().remove(r)

def renumber_rows(table):
    # Assumes we are renumbering a standard table
    tabname = str(table.name)
    idtag = sdmpy.sdm.decap(tabname) + 'Id'
    idmap = {}
    for (i,r) in enumerate(table):
        newid = tabname + '_%d'%i
        oldid = r.__getattr__(idtag)
        idmap[oldid] = newid
        r.__setattr__(idtag,newid)
    # Return the old->new mapping in case it's needed
    return idmap

def rename_rows(table,idmap,idtag='spectralWindowId'):
    for r in table:
        oldid = r.__getattr__(idtag)
        r.__setattr__(idtag,idmap[oldid])

# Not sure whether needed.  CASA seems to import OK without it, not sure
# if will be problems downstream.   Keeping original SpW ID names would
# presumably allow attaching original SysPower data if needed.  Looks like
# normal realfast SDMs keep the original SpW numbering so I think this setting
# (False) is correct:
redo_spws = False

# --- SpectralWindow table

clean_rows(sdm['SpectralWindow'],spwid_keep)
if redo_spws: spw_map = renumber_rows(sdm['SpectralWindow'])

# --- Feed table

clean_rows(sdm['Feed'],spwid_keep)
if redo_spws: rename_rows(sdm['Feed'],spw_map)

# --- Source table

clean_rows(sdm['Source'],spwid_keep)
if redo_spws: rename_rows(sdm['Source'],spw_map)

# --- DataDescription table

clean_rows(sdm['DataDescription'],spwid_keep)
if redo_spws: rename_rows(sdm['DataDescription'],spw_map)
dd_map = renumber_rows(sdm['DataDescription'])
# really this should translate not just assume ordering will be OK..
# but in practice this should be fine for these SDMs.
dd_keep = [str(r.dataDescriptionId) for r in sdm['DataDescription']]

# --- ConfigDescription table

# assume only one row..
cd = sdm['ConfigDescription'][0]
n = len(dd_keep)
cd.numDataDescription = n
cd.dataDescriptionId = ('1 %d ' % n) + ' '.join(dd_keep)
cd.switchCycleId = ('1 %d' % n) + ' SwitchCycle_0'*n

# --- Polarization table

# Assume this is the one to go, really could do this based on revised DD table..
r = sdm['Polarization']['Polarization_2']
r.getparent().remove(r)

print("Writing results to '%s'" % sdmout)
sdm.write(sdmout)
