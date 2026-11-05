import math

import pytest

from biomass_heat import thermo as th
from biomass_heat.agents import (Burner, Electrolyser, GasifierParams, Gasifier, GasStorage,
                                 methanation_equilibrium)
from biomass_heat.biomass import Biomass
from biomass_heat.configs import ALL_CONFIGS, natural_gas
from biomass_heat.plant import Plant
from biomass_heat.streams import GasStream


@pytest.mark.parametrize("key", list(ALL_CONFIGS))
def test_gasifier_element_and_energy_balance(key):
    cfg = ALL_CONFIGS[key]
    b = cfg.biomass
    gz = Gasifier(cfg.gasifier)
    feed = {"dry_kg_s": 0.85, "moisture_kg_s": 0.15}
    out = gz.process(feed, b)
    gas = out["gas"]
    a = {k: v * 0.85 for k, v in b.atoms_per_kg_dry().items()}
    f = 1 - out["recycle_frac"]
    at = gas.atoms()
    # carbon: gas + char (+ recycled share for allothermal mode)
    assert at["C"] / f + out["char_mol_s"] == pytest.approx(a["C"], rel=1e-9)
    steam = out["steam_kg_s"] / 18.015e-3
    h_in = a["H"] + 2 * (0.15 / 18.015e-3 + steam)
    assert at["H"] / f == pytest.approx(h_in, rel=1e-9)
    if cfg.gasifier.mode != "steam":
        o_in = a["O"] + 0.15 / 18.015e-3 + steam + 2 * out["o2_mol_s"]
        assert at["O"] == pytest.approx(o_in, rel=1e-9)


@pytest.mark.parametrize("key", ["air", "oxy_none"])
def test_gasifier_wgs_equilibrium_and_autothermal_balance(key):
    cfg = ALL_CONFIGS[key]
    out = Gasifier(cfg.gasifier).process({"dry_kg_s": 0.85, "moisture_kg_s": 0.15}, cfg.biomass)
    n = out["gas"].n
    K = th.equilibrium_constant({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, cfg.gasifier.T_out)
    assert n["CO2"] * n["H2"] / (n["CO"] * n["H2O"]) == pytest.approx(K, rel=1e-6)
    # energy: LHV in = LHV gas + sensible heat of gas + char LHV + loss (oxidant at 25 C)
    char_lhv = out["char_mol_s"] * -th.hf("CO2")
    sens = out["gas"].enthalpy - th.mix_h(n, th.T_REF)
    lhv_gas_wet = out["gas"].lhv_power
    steam_mol = out["steam_kg_s"] / 18.015e-3
    steam_sens = steam_mol * (th.h("H2O", cfg.gasifier.T_steam) - th.hf("H2O"))
    lhs = lhv_gas_wet + sens + char_lhv + out["q_loss_W"]
    assert lhs == pytest.approx(out["lhv_in_W"] + steam_sens, rel=2e-3)
    assert 0.2 < out["equivalence_ratio"] < 0.4


def test_air_gas_is_nitrogen_diluted_and_oxy_gas_is_not():
    air = Plant("air").chain(1.0)["cleaning"]["gas"].x(dry=True)
    oxy = Plant("oxy_none").chain(1.0)["cleaning"]["gas"].x(dry=True)
    assert air["N2"] > 0.4 and oxy.get("N2", 0) < 0.01
    assert Plant("oxy_none").chain(1.0)["cleaning"]["gas"].lhv_volumetric() > \
        1.5 * Plant("air").chain(1.0)["cleaning"]["gas"].lhv_volumetric()


def test_methanation_equilibrium_constants_satisfied():
    n0 = {"CO": 1.0, "CO2": 1.0, "H2": 7.0, "H2O": 0.1, "CH4": 0.2, "N2": 0.05}
    T, p = 573.15, 10e5
    eq = methanation_equilibrium(n0, T, p)
    tot = sum(eq.values())
    y = {k: v / tot for k, v in eq.items()}
    P = p / th.P_REF
    K1 = th.equilibrium_constant({"CO": -1, "H2": -3, "CH4": 1, "H2O": 1}, T)
    K2 = th.equilibrium_constant({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, T)
    assert y["CH4"] * y["H2O"] / (y["CO"] * y["H2"] ** 3) / P**2 == pytest.approx(K1, rel=1e-4)
    assert y["CO2"] * y["H2"] / (y["CO"] * y["H2O"]) == pytest.approx(K2, rel=1e-4)
    # element conservation
    for i in range(3):
        before = sum(v * th.ATOMS[k][i] for k, v in n0.items())
        after = sum(v * th.ATOMS[k][i] for k, v in eq.items())
        assert after == pytest.approx(before, rel=1e-8)


def test_methanation_pressure_and_temperature_trends():
    n0 = {"CO2": 1.0, "H2": 4.0}
    ch4 = lambda T, p: methanation_equilibrium(n0, T, p)["CH4"]
    assert ch4(573.15, 20e5) > ch4(573.15, 1e5)  # Le Chatelier: fewer moles
    assert ch4(573.15, 5e5) > ch4(873.15, 5e5)  # exothermic


def test_electrolyser_stoichiometry_and_limits():
    el = Electrolyser()
    out = el.process(55.0 * 3.6e6 / 3600)  # 55 kWh/h -> 1 kg/h H2
    assert out["h2_mol_s"] * 2.016e-3 * 3600 == pytest.approx(1.0, rel=1e-9)
    assert out["o2_mol_s"] == pytest.approx(out["h2_mol_s"] / 2)
    el2 = Electrolyser(type(el.p)(capacity_W=1e6, min_load=0.2))
    assert el2.process(5e6)["electricity_W"] == 1e6
    assert el2.process(0.1e6)["electricity_W"] == 0.0
    with pytest.raises(ValueError):
        el.process(-1.0)


def test_storage_limits():
    s = GasStorage("h2", capacity_mol=100.0, max_rate_mol_s=0.02)
    r = s.exchange(1.0, 3600)  # rate-limited
    assert r["flow_mol_s"] == pytest.approx(0.02)
    assert s.level == pytest.approx(72.0)
    r = s.exchange(0.02, 3600)  # capacity-limited
    assert s.level == pytest.approx(100.0)
    assert r["unserved_mol_s"] > 0
    r = s.exchange(-1.0, 3600)
    assert r["flow_mol_s"] == pytest.approx(-0.02)


def test_burner_methane_flame_temperature_and_available_heat():
    b = Burner(type(Burner().p)(excess_air=1.0))
    t_ad = b.adiabatic_flame_temperature(GasStream({"CH4": 1.0}, 298.15)) - 273.15
    # equilibrium value ~1950 C; complete combustion without dissociation is higher (~2050 C)
    assert 1950 < t_ad < 2100
    ng = natural_gas()
    burner = Burner()
    assert burner.available_heat_fraction(ng, 298.15 + 1e-9) == pytest.approx(1.0, abs=0.01)
    assert burner.available_heat_fraction(ng, 1273.15) < burner.available_heat_fraction(ng, 773.15)


def test_gasifier_param_validation():
    with pytest.raises(ValueError):
        GasifierParams(mode="plasma")
    with pytest.raises(ValueError):
        GasifierParams(ch4_carbon_fraction=0.99, carbon_conversion=0.9)
