import numpy as np
import pytest

from biomass_heat import CONFIGS, ALL_CONFIGS, GasStorage, Plant, PriceThreshold


@pytest.fixture(scope="module")
def kpis():
    return {k: Plant(k).design_point(10e6, 1273.15) for k in ALL_CONFIGS}


def test_design_point_meets_heat_duty(kpis):
    for k, r in kpis.items():
        assert r["available_heat_frac"] * r["fuel_lhv_MW"] == pytest.approx(10.0, rel=1e-6)


def test_fuel_quality_ranking(kpis):
    lhv = [kpis[k]["fuel_lhv_MJ_Nm3"] for k in CONFIGS]
    assert lhv == sorted(lhv)  # air < oxy < partial < full
    assert kpis["oxy_full"]["fuel_wobbe_MJ_Nm3"] > 45  # close to natural gas (~52.7)
    assert kpis["air"]["fuel_wobbe_MJ_Nm3"] < 10


def test_hydrogen_and_oxygen_logic(kpis):
    assert kpis["air"]["h2_kg_per_MWh"] == 0 and kpis["air"]["o2_kg_per_MWh"] == 0
    assert kpis["oxy_none"]["o2_asu_kg_h"] == pytest.approx(kpis["oxy_none"]["o2_demand_kg_h"])
    assert kpis["oxy_full"]["h2_kg_per_MWh"] > kpis["oxy_partial"]["h2_kg_per_MWh"] > 0
    # electrolyser O2 co-product exceeds the gasifier's demand when H2 is used for methanation
    for k in ("oxy_partial", "oxy_full"):
        assert kpis[k]["o2_coproduct_kg_h"] > kpis[k]["o2_demand_kg_h"]
        assert kpis[k]["o2_asu_kg_h"] == 0
    e = kpis["oxy_none_elec"]
    assert e["o2_coproduct_kg_h"] == pytest.approx(e["o2_demand_kg_h"], rel=1e-9)


def test_biomass_demand_falls_with_methanation(kpis):
    b = [kpis[k]["biomass_MWh_per_MWh"] for k in CONFIGS]
    assert b[0] > b[1] > b[2] > b[3]


def test_efficiencies_physical(kpis):
    for r in kpis.values():
        assert 0 < r["cold_gas_efficiency"] < 1
        assert 0 < r["fuel_conversion_eff"] < 1
        assert 0 < r["chain_eff_heat"] < 1


def test_simulation_matches_design_point():
    p = Plant("oxy_full")
    dp = p.design_point(10e6, 1273.15)
    df = p.simulate([10e6] * 24, [50.0] * 24)
    assert df["electricity_MW"].iloc[0] == pytest.approx(dp["electricity_MW"], rel=1e-6)
    assert df["h2_unserved_kg"].sum() == pytest.approx(0.0, abs=1e-9)
    assert df["electricity_cost_eur"].sum() == pytest.approx(dp["electricity_MW"] * 50 * 24, rel=1e-6)


def test_price_threshold_strategy_shifts_load_into_cheap_hours():
    p = Plant("oxy_full")
    dp = p.design_point(10e6, 1273.15)
    h2_mol_s = dp["h2_kg_h"] / 3600 / 2.016e-3
    prices = np.where(np.arange(48) % 24 < 12, 20.0, 120.0)
    store = GasStorage("h2_storage", capacity_mol=h2_mol_s * 3600 * 14, max_rate_mol_s=h2_mol_s * 3)
    cap = 2.2 * dp["el_electrolyser_MW"] * 1e6
    base = p.simulate([10e6] * 48, prices)
    smart = p.simulate([10e6] * 48, prices, strategy=PriceThreshold(20.0, cap), h2_storage=store)
    assert smart["electricity_cost_eur"].sum() < base["electricity_cost_eur"].sum()
    assert smart["h2_unserved_kg"].sum() == pytest.approx(0.0, abs=1e-6)


def test_infeasible_temperature_raises():
    with pytest.raises(ValueError):
        Plant("air").design_point(10e6, 1700 + 273.15)
