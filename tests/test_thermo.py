import math

import pytest

from biomass_heat import thermo as th
from biomass_heat.biomass import Biomass
from biomass_heat.streams import GasStream


@pytest.mark.parametrize("sp,lhv_kj", [("H2", 241.8), ("CO", 283.0), ("CH4", 802.3)])
def test_lhv_matches_standard_values(sp, lhv_kj):
    # standard heats of combustion (NIST/JANAF-consistent), +-1 kJ/mol
    assert th.lhv_molar(sp) / 1e3 == pytest.approx(lhv_kj, abs=1.0)


def test_formation_enthalpies():
    assert th.hf("CO2") / 1e3 == pytest.approx(-393.5, abs=0.2)
    assert th.hf("H2O") / 1e3 == pytest.approx(-241.8, abs=0.2)
    assert th.hf("H2") == pytest.approx(0.0, abs=1.0)


def test_wgs_equilibrium_constant_crosses_one_near_820C():
    # K_WGS = 1 at roughly 1090-1100 K (classic textbook value)
    k = lambda T: th.equilibrium_constant({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, T)
    assert k(1000) > 1 > k(1200)


def test_methane_properties_and_wobbe():
    ch4 = GasStream({"CH4": 1.0})
    assert ch4.lhv_volumetric() == pytest.approx(35.8, abs=0.2)  # MJ/Nm3
    assert ch4.hhv_volumetric() == pytest.approx(39.7, abs=0.2)
    assert ch4.wobbe_index() == pytest.approx(53.4, abs=0.4)


def test_solve_temperature_roundtrip():
    n = {"CO2": 1.0, "H2O": 2.0, "N2": 7.0}
    H = th.mix_h(n, 1500.0)
    assert th.solve_temperature(n, H) == pytest.approx(1500.0, abs=1e-3)


def test_p_sat_water():
    assert th.p_sat_water(373.15) == pytest.approx(101325, rel=0.01)
    assert th.p_sat_water(313.15) == pytest.approx(7384, rel=0.01)


def test_biomass_heating_values():
    b = Biomass()
    assert 19.0 < b.hhv_dry < 21.0  # typical wood, MJ/kg dry
    assert b.lhv_ar < b.lhv_dry < b.hhv_dry
    with pytest.raises(ValueError):
        Biomass(C=60)


def test_mixing_conserves_enthalpy():
    a = GasStream({"N2": 1.0}, 400.0)
    b = GasStream({"N2": 1.0}, 600.0)
    m = a + b
    assert m.enthalpy == pytest.approx(a.enthalpy + b.enthalpy, rel=1e-9)
    assert 495 < m.T < 505
