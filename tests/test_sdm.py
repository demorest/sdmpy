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

def test_tables(sdm):
    assert len(sdm.tables) == 34


def test_maintable(sdm):
    assert len(sdm['Main']) == 14
    assert sdm['Main'][0].numAntenna == 27


def test_bdfs(sdm):
    assert any([scan.bdf.exists for scan in sdm.scans()])


def test_meta(sdm):
    scan = sdm.scan(1)
    assert scan.numIntegration == 596
    assert len(scan.antennas) == 27
    assert len(scan.spws) == 8


def test_bdfmeta(sdm):
    bdf = list(sdm.scans())[6].bdf
    assert bdf.exists
    assert bdf.numIntegration == 1
    assert bdf.numAntenna == 27
    assert len(bdf.spws) == 8


def test_bdfdata0(sdm):
    bdf = list(sdm.scans())[6].bdf
    data = bdf.get_integration(0)
    assert data.data['crossData'].shape == (351, 512)


def test_bdfdata1(sdm):
    bdf = list(sdm.scans())[6].bdf
    data = bdf.get_data()
    assert data.shape == (1, 351, 8, 1, 32, 2)


def test_bdfdata2(sdm):
    bdf = list(sdm.scans())[6].bdf
    data = bdf.get_data(spwidx=0)
    assert data.shape == (1, 351, 1, 32, 2)


def test_bdfdata3(sdm):
    bdf = list(sdm.scans())[6].bdf
    data = bdf.get_data(fscrunch=True)
    assert data.shape == (1, 351, 8, 1, 2)


def test_spw(sdm):
    scan = sdm.scan(1)
    assert len(scan.spws) == 8
    spw = scan.spw(0)
    assert spw.chanWidth == 4000000.0
    assert str(spw.spectralWindowId) == 'SpectralWindow_0'
