"""Plant: scheduler and message bus that wires the agents of one process chain together.

Two entry points:

* :meth:`Plant.design_point` - steady-state mass and energy balance for a given process
  heat duty (all unit models are linear in throughput, so the chain is evaluated for a
  reference feed and scaled exactly).
* :meth:`Plant.simulate` - hourly time-stepping with a heat-demand profile, an electricity
  price series, H2/O2 storages and a pluggable electrolyser operating strategy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Optional, Protocol, Sequence

import pandas as pd

from . import thermo as th
from .agents import (MW, AirSeparationUnit, Biomass, BiomassFeed, Burner, Electrolyser,
                     GasCleaning, Gasifier, GasStorage, GridMarket, Methanation,
                     ProcessHeatDemand)
from .configs import ALL_CONFIGS, PlantConfig, natural_gas
from .streams import GasStream

KG_H = 3600.0  # kg/s -> kg/h


class Plant:
    def __init__(self, config: PlantConfig | str):
        if isinstance(config, str):
            config = ALL_CONFIGS[config]
        self.config = c = config
        self.feed = BiomassFeed(c.biomass)
        self.gasifier = Gasifier(c.gasifier)
        self.cleaning = GasCleaning(c.cleaning)
        self.methanation = Methanation(c.methanation, T_dry=c.cleaning.T_clean)
        self.electrolyser = Electrolyser(c.electrolyser)
        self.asu = AirSeparationUnit(c.asu)
        self.burner = Burner(c.burner)

    @property
    def agents(self):
        return [self.feed, self.gasifier, self.cleaning, self.methanation, self.electrolyser,
                self.asu, self.burner]

    # ------------------------------------------------------------------------------------
    def chain(self, feed_ar_kg_s: float = 1.0) -> Dict[str, Any]:
        """Steady-state flows through the chain for a given as-received biomass feed."""
        c = self.config
        feed = self.feed.process(feed_ar_kg_s)
        gz = self.gasifier.process(feed, c.biomass)
        cl = self.cleaning.process(gz["gas"])
        h2_meth = self.methanation.h2_demand(cl["gas"])
        me = self.methanation.process(cl["gas"], h2_meth)
        fuel = me["gas"]
        o2_need = gz["o2_mol_s"] if c.gasifier.mode == "oxyfuel" else 0.0

        h2_el = h2_meth
        h2_blend = 0.0
        if c.o2_supply == "electrolyser" and c.blend_surplus_h2:
            h2_el = max(h2_meth, 2 * o2_need)
            h2_blend = h2_el - h2_meth
        el = self.electrolyser.process(self.electrolyser.power_for_h2(h2_el)) if h2_el > 0 else \
            self.electrolyser.process(0.0)
        if h2_blend > 0:
            fuel = fuel + GasStream({"H2": h2_blend}, fuel.T, fuel.p)

        o2_from_el = min(el["o2_mol_s"], o2_need) if c.o2_supply == "electrolyser" else 0.0
        o2_asu = o2_need - o2_from_el if c.gasifier.mode == "oxyfuel" else 0.0
        asu = self.asu.process(o2_asu)
        o2_surplus = el["o2_mol_s"] - o2_from_el

        steam_mol = gz["steam_kg_s"] * 1000 / th.MOLAR_MASS["H2O"]
        steam_heat = steam_mol * (th.h("H2O", c.gasifier.T_steam) - (th.hf("H2O") - th.L_H2O_25C))
        return {"feed": feed, "gasifier": gz, "cleaning": cl, "methanation": me,
                "electrolyser": el, "asu": asu, "fuel": fuel, "o2_need_mol_s": o2_need,
                "o2_from_electrolyser_mol_s": o2_from_el, "o2_asu_mol_s": o2_asu,
                "o2_surplus_mol_s": o2_surplus, "h2_methanation_mol_s": h2_meth,
                "h2_blend_mol_s": h2_blend, "steam_heat_W": steam_heat}

    # ------------------------------------------------------------------------------------
    def design_point(self, heat_W: float = 10 * MW, T_process: float = 1000 + 273.15,
                     flue_exit_approach_K: float = 0.0) -> Dict[str, Any]:
        """Size the chain for a process heat duty and return key performance indicators."""
        T_exit = T_process + flue_exit_approach_K
        ref = self.chain(1.0)
        q_ref = self.burner.available_heat(ref["fuel"], T_exit)
        if q_ref <= 0:
            raise ValueError("fuel gas cannot supply heat at this process temperature "
                             "(adiabatic flame temperature below flue-gas exit temperature)")
        s = heat_W / q_ref
        r = self.chain(s)
        return self._kpis(r, heat_W, T_exit)

    def _kpis(self, r: Dict[str, Any], heat_W: float, T_exit: float) -> Dict[str, Any]:
        c = self.config
        gz, cl, me, el, asu, fuel = (r["gasifier"], r["cleaning"], r["methanation"],
                                     r["electrolyser"], r["asu"], r["fuel"])
        syngas = cl["gas"]
        bio_W = r["feed"]["lhv_W"]
        h2_lhv_W = el["h2_mol_s"] * th.lhv_molar("H2")
        elec = {"electrolyser": el["electricity_W"], "asu": asu["electricity_W"],
                "compression": me["electricity_W"], "blower": cl["electricity_W"]}
        e_tot = sum(elec.values())
        mwh = heat_W  # per-MWh normalisation: X_W / heat_W = MWh_X per MWh_heat
        to_kg_h = lambda mol_s, sp: mol_s * th.MOLAR_MASS[sp] / 1000 * KG_H
        ng = natural_gas()
        ng_frac = self.burner.available_heat_fraction(ng, T_exit)
        ng_lhv_per_heat = 1 / ng_frac
        ng_co2_t_per_mwh_heat = (self.burner.co2_kg_s(ng) / ng.lhv_power) * ng_lhv_per_heat * 3.6e6
        heat_rec = cl["heat_recovered_W"] + me["heat_recovered_W"]
        out = {
            "config": c.key, "label": c.label,
            "process_heat_MW": heat_W / MW, "flue_exit_C": T_exit - 273.15,
            "biomass_feed_t_h_ar": r["feed"]["dry_kg_s"] / (1 - c.biomass.moisture / 100) * 3.6,
            "biomass_lhv_MW": bio_W / MW,
            "equivalence_ratio": gz.get("equivalence_ratio", 0.0),
            "cold_gas_efficiency": syngas.lhv_power / bio_W,
            "syngas_lhv_MJ_Nm3": syngas.lhv_volumetric(),
            "fuel_lhv_MJ_Nm3": fuel.lhv_volumetric(),
            "fuel_hhv_MJ_Nm3": fuel.hhv_volumetric(),
            "fuel_wobbe_MJ_Nm3": fuel.wobbe_index(),
            "fuel_rel_density": fuel.relative_density(),
            "fuel_lhv_MW": fuel.lhv_power / MW,
            "adiabatic_flame_C": self.burner.adiabatic_flame_temperature(fuel) - 273.15,
            "available_heat_frac": heat_W / fuel.lhv_power,
            "h2_kg_h": to_kg_h(el["h2_mol_s"], "H2"),
            "o2_demand_kg_h": to_kg_h(r["o2_need_mol_s"], "O2"),
            "o2_coproduct_kg_h": to_kg_h(el["o2_mol_s"], "O2"),
            "o2_surplus_kg_h": to_kg_h(r["o2_surplus_mol_s"], "O2"),
            "o2_asu_kg_h": to_kg_h(r["o2_asu_mol_s"], "O2"),
            **{f"el_{k}_MW": v / MW for k, v in elec.items()},
            "electricity_MW": e_tot / MW,
            "heat_recovered_MW": heat_rec / MW,
            "steam_raising_MW": r["steam_heat_W"] / MW,
            # normalised per MWh of process heat
            "biomass_MWh_per_MWh": bio_W / mwh,
            "electricity_MWh_per_MWh": e_tot / mwh,
            "h2_kg_per_MWh": to_kg_h(el["h2_mol_s"], "H2") / (heat_W / MW),
            "o2_kg_per_MWh": to_kg_h(r["o2_need_mol_s"], "O2") / (heat_W / MW),
            "fuel_conversion_eff": fuel.lhv_power / (bio_W + h2_lhv_W),
            "chain_eff_heat": heat_W / (bio_W + e_tot),
            "ng_available_heat_frac": ng_frac,
            "ng_MWh_per_MWh": ng_lhv_per_heat,
            "ng_avoided_co2_t_per_MWh": ng_co2_t_per_mwh_heat,
        }
        for sp in ("H2", "CO", "CO2", "CH4", "N2", "Ar"):
            out[f"syngas_{sp}_vol_dry"] = syngas.x(dry=True).get(sp, 0.0) * 100
            out[f"fuel_{sp}_vol_dry"] = fuel.x(dry=True).get(sp, 0.0) * 100
        out["syngas_H2O_vol_wet"] = syngas.x().get("H2O", 0.0) * 100
        out["fuel_H2O_vol_wet"] = fuel.x().get("H2O", 0.0) * 100
        return out

    # ------------------------------------------------------------------------------------
    def simulate(self, demand_W: Sequence[float], prices_eur_mwh: Optional[Sequence[float]] = None,
                 T_process: float = 1000 + 273.15, strategy: Optional["Strategy"] = None,
                 h2_storage: Optional[GasStorage] = None, o2_storage: Optional[GasStorage] = None,
                 dt: float = 3600.0) -> pd.DataFrame:
        """Hourly agent-based simulation.

        Each step: the demand agent posts the heat duty, the chain is sized for it, the
        strategy sets the electrolyser power, storages absorb the H2/O2 mismatch, unserved
        O2 is bought from the ASU, surplus O2 is vented and the grid agent books the
        electricity cost. A H2 shortfall (empty storage) is reported as ``h2_unserved_kg``.
        """
        strategy = strategy or FollowDemand()
        demand = ProcessHeatDemand(demand_W, T_process)
        grid = GridMarket(prices_eur_mwh if prices_eur_mwh is not None else [0.0])
        h2s = h2_storage or GasStorage("h2_storage")
        o2s = o2_storage or GasStorage("o2_storage")
        for a in [*self.agents, demand, grid, h2s, o2s]:
            a.reset()
        uses_o2_el = self.config.o2_supply == "electrolyser"
        ref = self.chain(1.0)
        q_ref = self.burner.available_heat(ref["fuel"], demand.T_flue_exit)
        rows = []
        for t in range(len(demand.profile)):
            d = demand.step(t, dt, {})
            s = d["demand_W"] / q_ref
            r = self.chain(s) if s > 0 else None
            h2_need = (r["h2_methanation_mol_s"] + r["h2_blend_mol_s"]) if r else 0.0
            o2_need = r["o2_need_mol_s"] if (r and uses_o2_el) else 0.0
            price = grid.price(t)
            p_el = strategy.electrolyser_power(t, self, h2_need, o2_need, price, h2s, o2s)
            el = self.electrolyser.step(t, dt, {"power_W": p_el})
            h2x = h2s.step(t, dt, {"net_in_mol_s": el["h2_mol_s"] - h2_need})
            o2_net = (el["o2_mol_s"] if uses_o2_el else 0.0) - o2_need
            o2x = o2s.step(t, dt, {"net_in_mol_s": o2_net})
            o2_short = max(0.0, -o2x["unserved_mol_s"])
            o2_vent = max(0.0, o2x["unserved_mol_s"])
            o2_asu_base = r["o2_asu_mol_s"] if (r and not uses_o2_el) else 0.0
            asu = self.asu.step(t, dt, {"o2_mol_s": o2_asu_base + o2_short})
            other_el = (r["methanation"]["electricity_W"] + r["cleaning"]["electricity_W"]) if r else 0.0
            e_tot = el["electricity_W"] + asu["electricity_W"] + other_el \
                + h2x["electricity_W"] + o2x["electricity_W"]
            grid.step(t, dt, {"electricity_W": e_tot})
            rows.append({
                "t": t, "price_eur_mwh": price, "heat_MW": d["demand_W"] / MW,
                "biomass_MW": (r["feed"]["lhv_W"] if r else 0.0) / MW,
                "electrolyser_MW": el["electricity_W"] / MW, "asu_MW": asu["electricity_W"] / MW,
                "electricity_MW": e_tot / MW,
                "h2_need_kg": h2_need * 2.016e-3 * dt,
                "h2_produced_kg": el["h2_mol_s"] * 2.016e-3 * dt,
                "h2_unserved_kg": max(0.0, -h2x["unserved_mol_s"]) * 2.016e-3 * dt,
                "h2_spill_kg": max(0.0, h2x["unserved_mol_s"]) * 2.016e-3 * dt,
                "h2_storage_kg": h2s.level * 2.016e-3,
                "o2_storage_kg": o2s.level * 31.998e-3,
                "o2_vent_kg": o2_vent * 31.998e-3 * dt,
                "electricity_cost_eur": e_tot / MW * price * dt / 3600,
            })
        return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------
class Strategy(Protocol):
    def electrolyser_power(self, t: int, plant: Plant, h2_need_mol_s: float,
                           o2_need_mol_s: float, price: float, h2_storage: GasStorage,
                           o2_storage: GasStorage) -> float: ...


@dataclass
class FollowDemand:
    """Baseline: run the electrolyser exactly at the hourly H2 (or O2) requirement."""

    def electrolyser_power(self, t, plant, h2_need_mol_s, o2_need_mol_s, price, h2_storage,
                           o2_storage):
        h2 = max(h2_need_mol_s, 2 * o2_need_mol_s if plant.config.blend_surplus_h2 else 0.0)
        return plant.electrolyser.power_for_h2(h2)


@dataclass
class PriceThreshold:
    """Simple illustration of price-dependent operation: run at full capacity when the price
    is below ``threshold`` (filling the H2 store), otherwise draw from storage and only run
    for the part storage cannot cover. Optimised strategies follow in the dispatch repo."""

    threshold_eur_mwh: float
    capacity_W: float

    def electrolyser_power(self, t, plant, h2_need_mol_s, o2_need_mol_s, price, h2_storage,
                           o2_storage):
        need_W = plant.electrolyser.power_for_h2(h2_need_mol_s)
        if price <= self.threshold_eur_mwh:
            return max(self.capacity_W, need_W)
        from_store = min(h2_storage.level / 3600.0, h2_need_mol_s, h2_storage.max_rate)
        return plant.electrolyser.power_for_h2(h2_need_mol_s - from_store)
