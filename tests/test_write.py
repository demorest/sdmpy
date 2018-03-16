import sdmpy
import pytest
import os.path
import numpy as np

_install_dir = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture(scope="module")
def sdm(request):
    sdmfile = os.path.join(_install_dir,
                           'data/16A-459_TEST_1hr_000.57633.66130137732.scan7.cut1')

    return sdmpy.SDM(sdmfile)


def test_read(sdm):
    bdf = list(sdm.scans())[6].bdf
    data = bdf.get_data()
    nint, nbl, nspw, numBin, nchan, npol = data.shape

    spws = [sdmpy.bdf.BDFSpectralWindow(None, numBin=numBin,
                                        numSpectralPoint=nchan,
                                        sw=1, swbb='AC_8BIT', npol=2)]

    uid = 'uid:///evla/test/0123456789'
    w = sdmpy.bdf.BDFWriter('.', start_mjd=0, uid=uid,
                            num_antenna=27, spws=spws, scan_idx=1,
                            corr_mode='c')

    dat = {}
    w.write_header()
    for i in range(nint):
        dat['crossData'] = data[i]
        ts = 0
        w.write_integration(mjd=ts, interval=1, data=dat)
    w.close()
