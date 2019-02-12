[![Build Status](https://travis-ci.org/realfastvla/sdmpy.svg?branch=master)](https://travis-ci.org/realfastvla/sdmpy)
[![codecov](https://codecov.io/gh/realfastvla/sdmpy/branch/master/graph/badge.svg)](https://codecov.io/gh/realfastvla/sdmpy)

sdmpy -- P. Demorest, 2015/03

Simple classes for reading VLA/ALMA Science Data Model (SDM) XML and
data files.  The structure of the SDM is a set of tables, which contain
rows of data.  The main goal of this is to have a straightforward,
relatively dependency-free way of doing simple or low-level stuff with
SDMs (e.g., for situations where invoking CASA is unnecessary).

Most of the sdmpy functionality is accessed via one of the following
three classes:

  1. SDM class -- for reading or editing SDM XML tables.

  2. BDF class -- for reading data from Binary Data Format (BDF) files.
     This code can be used directly on standalone BDF files, there is no
     need for a full SDM to be available.

  3. Scan class -- provides higher-level access to a full SDM+BDF data
     set, including most commonly needed metadata (times, frequencies,
     source names, antenna names, etc.)

There are also a handful of command-line programs making use of sdmpy in
the scripts/ directory.  Most can be run with '-h' option for usage
information.

---------------------------------------------------------------------
Example ipython session illustrating use of the SDM class:

```
In [1]: import sdmpy

In [2]: s = sdmpy.SDM('15A-105_sb30463808_1.57082.53091797454')

In [3]: s.tables
Out[3]:
['SpectralWindow',
 'CalData',
 'Antenna',
 'Subscan',
 'SwitchCycle',
 'Polarization',
 'Source',
 'State',
 'Station',
 'Main',
 'Flag',
 'ExecBlock',
 'SBSummary',
 'ConfigDescription',
 'Receiver',
 'Processor',
 'CorrelatorMode',
 'Feed',
 'PointingModel',
 'Field',
 'Scan',
 'CalDevice',
 'Weather',
 'CalReduction',
 'DataDescription']

In [4]: s['Field'][0].keys
Out[4]:
['code',
 'referenceDir',
 'phaseDir',
 'sourceId',
 'delayDir',
 'fieldName',
 'time',
 'fieldId',
 'numPoly']

In [5]: s['Field'][0].fieldName
Out[5]: 'J1728+0427'

In [6]: for row in s['Field']: print row.fieldName
J1728+0427
J1713+0747
J1924-2914
J1909-3744
1411+522=3C295

In [7]: len(s['Main'])
Out[7]: 24
```

---------------------------------------------------------------------
ipython session showing use of the BDF class to access visibility data.
Note that this has been developed mainly for VLA/WIDAR data and is
likely not complete for use with ALMA data.  Example usage (see
docstrings for more info):

```
In [1]: import sdmpy

In [2]: b = sdmpy.BDF('uid____evla_bdf_1433189755525')

In [3]: b.numAntenna
Out[3]: 25

In [4]: b.numBaseline
Out[4]: 300

In [5]: b.numIntegration
Out[5]: 180

In [6]: b.basebands
Out[6]: ['AC_8BIT', 'BD_8BIT']

In [7]: i = b.get_integration(17)    # Get integration number 17

In [8]: d = i.get_data(0)            # Get spw 0 data

In [9]: d.shape
Out[9]: (300, 1, 64, 4)

In [10]: d[10,0,:,0]
Out[10]:
array([ -1.92527496e-03+0.00236825j,  -1.01333787e-03-0.0019273j ,
         3.88602726e-04-0.00326904j,   3.37082194e-04+0.0035117j ,
         4.74026240e-03-0.00303986j,  -1.39878283e-03+0.00315546j,
         4.75925812e-03+0.00384983j,  -6.46819070e-04-0.00230185j,
         1.04035577e-03-0.0033494j ,   4.77438327e-04+0.00407107j,
        -7.42492499e-03-0.00551728j,  -4.23104153e-04-0.00136438j,
        -1.09272427e-04-0.00189404j,   2.70528137e-04+0.00265754j,
         9.03953623e-04+0.00172536j,  -2.13428610e-03+0.00210991j,
         4.50677879e-04+0.00092774j,   2.75413692e-03+0.00105019j,
         2.04323884e-03+0.00331562j,  -4.14338522e-03+0.00393523j,
         2.90520466e-03-0.00068069j,  -7.84504786e-03+0.0031926j ,
        -1.59289513e-03-0.00068654j,  -8.57894833e-04+0.00268888j,
        -4.27102018e-03+0.00065515j,   3.69449146e-03-0.00241908j,
         4.95729037e-04+0.00161663j,  -1.60046294e-03-0.00291464j,
        -2.17019930e-03-0.00133767j,  -3.01234191e-04-0.00817809j,
         5.38951717e-05-0.00192556j,  -5.56168333e-03-0.00130252j,
         5.06524229e-05+0.00059094j,   1.00876275e-03+0.00976392j,
        -6.54946594e-03+0.00503897j,  -2.47276877e-03-0.0019854j ,
         3.56894545e-03+0.00648416j,  -1.66386040e-03-0.00040906j,
         3.25945253e-03-0.0011268j ,   6.37905393e-03+0.00060488j,
         2.66728364e-03-0.00056195j,   2.84236716e-03+0.00086546j,
         3.23694502e-03+0.00111199j,   6.35888521e-03+0.00346222j,
        -3.01071745e-03+0.0015606j ,   6.97700772e-04+0.00021654j,
        -5.21393435e-04+0.00208696j,   2.45397002e-03-0.00393461j,
        -1.13227195e-03+0.00167564j,  -2.68122321e-03-0.00271649j,
        -1.96605665e-03-0.00489836j,   2.49449024e-03+0.0007483j ,
        -7.28207640e-04+0.00322795j,  -3.14944983e-03-0.00022846j,
         9.21148807e-04+0.00533055j,  -9.83621459e-04+0.00440285j,
        -4.31790622e-03-0.00243727j,  -4.21724096e-03-0.00223194j,
         1.09594746e-03+0.00397917j,  -2.94885552e-03-0.00085458j,
         7.17548071e-04-0.00044755j,  -2.99423118e-05-0.00164009j,
         1.18584884e-03+0.00028139j,  -5.73034631e-04+0.00078536j], dtype=complex64)
```

---------------------------------------------------------------------
ipython session showing use of the Scan class:

```
In [1]: import sdmpy

In [2]: sdm = sdmpy.SDM('16B-248_TEST_L.57686.804786203706')

In [3]: for scan in sdm.scans(): print scan.idx,scan.source,scan.intents
1 1331+305=3C286 ['SYSTEM_CONFIGURATION']
2 1331+305=3C286 ['CALIBRATE_DELAY', 'CALIBRATE_FLUX',
'CALIBRATE_POL_ANGLE', 'CALIBRATE_BANDPASS']
3 J1626-2951 ['CALIBRATE_AMPLI', 'CALIBRATE_PHASE']
4 J1600-3053 ['OBSERVE_TARGET']
5 J1626-2951 ['CALIBRATE_AMPLI', 'CALIBRATE_PHASE']
6 J1600-3053 ['OBSERVE_TARGET']
7 J1626-2951 ['CALIBRATE_AMPLI', 'CALIBRATE_PHASE']

In [4]: s = sdm.scan(3)

In [5]: s.freqs()
Out[5]:
array([[  9.88000000e+08,   9.89000000e+08,   9.90000000e+08, ...,
          1.11300000e+09,   1.11400000e+09,   1.11500000e+09],
       [  1.11600000e+09,   1.11700000e+09,   1.11800000e+09, ...,
          1.24100000e+09,   1.24200000e+09,   1.24300000e+09],
       [  1.24400000e+09,   1.24500000e+09,   1.24600000e+09, ...,
          1.36900000e+09,   1.37000000e+09,   1.37100000e+09],
       ...,
       [  1.62800000e+09,   1.62900000e+09,   1.63000000e+09, ...,
          1.75300000e+09,   1.75400000e+09,   1.75500000e+09],
       [  1.75600000e+09,   1.75700000e+09,   1.75800000e+09, ...,
          1.88100000e+09,   1.88200000e+09,   1.88300000e+09],
       [  1.88400000e+09,   1.88500000e+09,   1.88600000e+09, ...,
          2.00900000e+09,   2.01000000e+09,   2.01100000e+09]])

In [6]: s.startMJD
Out[6]: 57686.808208912036

In [7]: b = s.bdf        # Get BDF object for this scan

In [8]: d = b.get_data() # Read full visibility data array for the scan
```

---------------------------------------------------------------------
ipython session showing how to read the SysPower table:

```
In [1]: import sdmpy

In [2]: import sdmpy.bintab

In [3]: s = sdmpy.SDM('17B-390_TEST_004.57994.73380716435')

In [4]: sp = sdmpy.bintab.unpacker(s['SysPower'])

In [5]: sp.unpack()

In [6]: sp.row
Out[6]:
rec.array([ ('Antenna_12', 'SpectralWindow_0', 0, 5010745004500000768, 1000000512, 2, [-0.95052397, -1.26523304], [  8.90080357,  11.66342068], [ 0.015625,  0.015625]),
 ('Antenna_12', 'SpectralWindow_1', 0, 5010745004500000768, 1000000512, 2, [ 1.625651  ,  1.82755697], [ 30.64261436,  27.64795876], [ 0.015625,  0.015625]),
 ('Antenna_12', 'SpectralWindow_2', 0, 5010745004500000768, 1000000512, 2, [-0.40402299, -0.38874799], [ 17.41805077,  14.52512169], [ 0.015625,  0.015625]),
 ...,
 ('Antenna_7', 'SpectralWindow_5', 0, 5010752180500000768, 1000000512, 2, [ 0.26653701,  0.311874  ], [ 21.20673752,  24.39814377], [ 0.022736,  0.022736]),
 ('Antenna_7', 'SpectralWindow_6', 0, 5010752180500000768, 1000000512, 2, [ 0.270816  ,  0.28408101], [ 20.47394753,  22.51699257], [ 0.020782,  0.022675]),
 ('Antenna_7', 'SpectralWindow_7', 0, 5010752180500000768, 1000000512, 2, [ 0.31178999,  0.27066699], [ 22.78490829,  22.75772285], [ 0.021698,  0.021759])],
          dtype=[('antennaId', 'S32'), ('spectralWindowId', 'S32'), ('feedId', '<i4'), ('timeMid', '<i8'), ('interval', '<i8'), ('numReceptor', '<i4'), ('switchedPowerDifference', '<f4', (2,)), ('switchedPowerSum', '<f4', (2,)), ('requantizerGain', '<f4', (2,))])
```
