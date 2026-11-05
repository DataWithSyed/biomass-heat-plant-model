"""Agents of the biomass process-heat plant.

Each agent is a self-contained unit with its own parameters and state. It exposes

* ``process(...)`` - the steady-state physics (mass and energy balance) of the unit, and
* ``step(t, dt, msg)`` - the time-step behaviour used by :class:`biomass_heat.plant.Plant`,
  which acts as the scheduler and message bus between agents.

Units: mol/s for gas flows, kg/s for solids/water, W for power, K for temperature, Pa for
pressure, unless the name says otherwise.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

from . import thermo as th
from .biomass import Biomass
from .streams import GasStream
from .units import compression_power, condense_water, sensible_heat

MW = 1e6
KWH = 3.6e6  # J


class Agent:
    """Base class: name, a key-value log of the last step, and a cumulative ledger."""

    name: str = "agent"

    def __init__(self, name: Optional[str] = None) -> None:
        if name:
            self.name = name
        self.last: Dict[str, float] = {}
        self.ledger: Dict[str, float] = {}

    def book(self, dt: float, **flows: float) -> None:
        """Record per-step values and add flow*dt to the cumulative ledger."""
        self.last = dict(flows)
        for k, v in flows.items():
            self.ledger[k] = self.ledger.get(k, 0.0) + v * dt

    def reset(self) -> None:
        self.last, self.ledger = {}, {}

    def step(self, t: int, dt: float, msg: Dict[str, Any]) -> Dict[str, Any]:  # pragma: no cover
        raise NotImplementedError


# ------------------------------------------------------------------------------------------
class BiomassFeed(Agent):
    """Supplies as-received biomass (dry matter + moisture) at the rate requested."""

    name = "biomass_feed"

    def __init__(self, biomass: Biomass = Biomass(), name: Optional[str] = None) -> None:
        super().__init__(name)
        self.biomass = biomass

    def process(self, feed_ar_kg_s: float) -> Dict[str, float]:
        w = self.biomass.moisture / 100
        return {"dry_kg_s": feed_ar_kg_s * (1 - w), "moisture_kg_s": feed_ar_kg_s * w,
                "lhv_W": feed_ar_kg_s * self.biomass.lhv_ar * MW}

    def step(self, t, dt, msg):
        out = self.process(msg["feed_ar_kg_s"])
        self.book(dt, feed_ar_kg_s=msg["feed_ar_kg_s"], biomass_lhv_W=out["lhv_W"])
        return out


# ------------------------------------------------------------------------------------------
@dataclass
class GasifierParams:
    """Engineering assumptions of the gasifier (see docs/assumptions.md)."""

    mode: str = "air"  # "air" | "oxyfuel" | "steam" (allothermal, dual fluidised bed)
    T_out: float = 850 + 273.15  # gasification / outlet temperature
    p: float = th.P_REF
    carbon_conversion: float = 0.95  # fraction of biomass C converted to gas
    ch4_carbon_fraction: float = 0.08  # fraction of biomass C leaving as CH4 (non-equilibrium)
    heat_loss_frac: float = 0.02  # heat loss as fraction of biomass LHV input
    steam_to_dry_biomass: float = 0.0  # kg steam / kg dry biomass
    T_steam: float = 300 + 273.15
    T_oxidant: float = th.T_REF
    o2_purity: float = 0.995  # oxyfuel only, balance Ar (electrolysis O2 is ~pure)
    # allothermal (steam) mode: heat from burning char (+ product gas if needed) in a combustor
    combustor_heat_transfer_eff: float = 0.80

    def __post_init__(self) -> None:
        if self.mode not in ("air", "oxyfuel", "steam"):
            raise ValueError("mode must be 'air', 'oxyfuel' or 'steam'")
        if not 0 < self.carbon_conversion <= 1:
            raise ValueError("carbon_conversion must be in (0, 1]")
        if not 0 <= self.ch4_carbon_fraction < self.carbon_conversion:
            raise ValueError("ch4_carbon_fraction must be < carbon_conversion")


class Gasifier(Agent):
    """Quasi-equilibrium gasifier.

    * Elemental C, H, O, N balances; unconverted carbon leaves as char.
    * CH4 is fixed as a fraction of the biomass carbon (pure equilibrium at 800-900 degC
      predicts almost no CH4, while real fluidised beds produce several vol-%; this
      non-equilibrium correction is common in equilibrium models, see Puig-Arnavat et al. 2010).
    * CO/CO2/H2/H2O are set by the water-gas-shift equilibrium at the outlet temperature.
    * Autothermal modes (air, oxyfuel): the oxidant flow is solved so that the energy balance
      closes at T_out (incl. heat loss).
    * Allothermal steam mode: no oxidant; the heat deficit is covered by burning the char
      and, if insufficient, a fraction of the product gas in a separate combustor.
    """

    name = "gasifier"

    def __init__(self, params: GasifierParams = GasifierParams(), name: Optional[str] = None):
        super().__init__(name)
        self.p = params

    # -- composition for a given oxygen supply ---------------------------------------------
    def _composition(self, el: Dict[str, float], n_ch4: float, n_inert: Dict[str, float]):
        T = self.p.T_out
        K = th.equilibrium_constant({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, T)
        A = el["C"] - n_ch4  # CO + CO2
        B = (el["H"] - 4 * n_ch4) / 2  # H2 + H2O
        O = el["O"]  # CO + 2 CO2 + H2O
        if A < 0 or B < 0 or O - A < 0:
            raise ValueError("infeasible elemental balance (not enough O or H)")
        # (1-K) y^2 + (B - O + A + K O) y - K A (O - A) = 0, y = CO2
        a, b, c = 1 - K, B - O + A + K * O, -K * A * (O - A)
        lo, hi = max(0.0, O - A - B), min(A, O - A)
        if abs(a) < 1e-12:
            y = -c / b
        else:
            disc = math.sqrt(max(b * b - 4 * a * c, 0.0))
            roots = [(-b + disc) / (2 * a), (-b - disc) / (2 * a)]
            cand = [r for r in roots if lo - 1e-9 <= r <= hi + 1e-9]
            if not cand:
                raise ValueError("no physical WGS root")
            y = cand[0]
        y = min(max(y, lo), hi)
        n = {"CO": A - y, "CO2": y, "H2O": O - A - y, "H2": B - (O - A - y), "CH4": n_ch4}
        n.update({k: v for k, v in n_inert.items() if v > 0})
        return GasStream(n, T, self.p.p)

    def process(self, feed: Dict[str, float], biomass: Biomass) -> Dict[str, Any]:
        p = self.p
        dry = feed["dry_kg_s"]
        a = {k: v * dry for k, v in biomass.atoms_per_kg_dry().items()}
        h2o_moist = feed["moisture_kg_s"] / th.MOLAR_MASS["H2O"] * 1000
        h2o_steam = p.steam_to_dry_biomass * dry / th.MOLAR_MASS["H2O"] * 1000
        c_gas = a["C"] * p.carbon_conversion
        n_char = a["C"] - c_gas
        n_ch4 = a["C"] * p.ch4_carbon_fraction
        lhv_in = dry * biomass.lhv_dry * MW - h2o_moist * th.L_H2O_25C
        q_loss = p.heat_loss_frac * lhv_in
        h_liq = th.hf("H2O") - th.L_H2O_25C
        H_in0 = (dry * biomass.formation_enthalpy_dry() + h2o_moist * h_liq
                 + h2o_steam * th.h("H2O", p.T_steam))

        def build(n_o2: float):
            if p.mode == "air":
                n_n2 = n_o2 * th.AIR["N2"] / th.AIR["O2"]
                n_ar = n_o2 * th.AIR["Ar"] / th.AIR["O2"]
            else:
                n_n2, n_ar = 0.0, n_o2 * (1 - p.o2_purity) / p.o2_purity
            el = {"C": c_gas, "H": a["H"] + 2 * (h2o_moist + h2o_steam),
                  "O": a["O"] + h2o_moist + h2o_steam + 2 * n_o2}
            inert = {"N2": a["N"] / 2 + n_n2, "Ar": n_ar}
            gas = self._composition(el, n_ch4, inert)
            H_ox = n_o2 * th.h("O2", p.T_oxidant) + n_n2 * th.h("N2", p.T_oxidant) \
                + n_ar * th.h("Ar", p.T_oxidant)
            resid = H_in0 + H_ox - gas.enthalpy - q_loss  # char at 25 degC: h = 0
            return gas, resid, n_n2, n_ar

        info: Dict[str, Any] = {"lhv_in_W": lhv_in, "q_loss_W": q_loss, "char_mol_s": n_char,
                                "steam_kg_s": h2o_steam * th.MOLAR_MASS["H2O"] / 1000,
                                "o2_mol_s": 0.0, "air_N2_mol_s": 0.0, "recycle_frac": 0.0,
                                "combustor_heat_W": 0.0}
        if p.mode in ("air", "oxyfuel"):
            n_o2_stoich = a["C"] + a["H"] / 4 - a["O"] / 2
            lo, hi = 0.0, n_o2_stoich
            if build(lo)[1] > 0:
                raise ValueError("gasifier needs no oxidant at these settings (check inputs)")
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if build(mid)[1] < 0:
                    lo = mid
                else:
                    hi = mid
            n_o2 = 0.5 * (lo + hi)
            gas, resid, n_n2, n_ar = build(n_o2)
            info.update(o2_mol_s=n_o2, air_N2_mol_s=n_n2, ar_mol_s=n_ar,
                        equivalence_ratio=n_o2 / n_o2_stoich)
        else:  # allothermal steam gasification (dual fluidised bed)
            gas, resid, _, _ = build(0.0)
            deficit = -resid  # W that must be supplied by the combustor
            char_heat = n_char * (-th.hf("CO2")) * p.combustor_heat_transfer_eff
            extra = max(0.0, deficit - char_heat)
            # fraction f of the product gas burnt in the combustor (LHV basis)
            f = min(1.0, extra / (gas.lhv_power * p.combustor_heat_transfer_eff))
            gas = gas.scaled(1 - f)
            info.update(equivalence_ratio=0.0, recycle_frac=f, combustor_heat_W=deficit,
                        char_heat_available_W=char_heat)
        info["cold_gas_efficiency"] = gas.lhv_power / lhv_in
        return {"gas": gas, **info}

    def step(self, t, dt, msg):
        out = self.process(msg["feed"], msg["biomass"])
        self.book(dt, syngas_lhv_W=out["gas"].lhv_power, o2_mol_s=out["o2_mol_s"])
        return out


# ------------------------------------------------------------------------------------------
@dataclass
class CleaningParams:
    T_clean: float = 40 + 273.15  # gas temperature after cooling/scrubbing
    T_heat_recovery_min: float = 200 + 273.15  # syngas cooler outlet (above tar dew point)
    heat_recovery_eff: float = 0.85
    pressure_drop: float = 0.2e5  # Pa, recompressed by a blower
    blower_eta: float = 0.70


class GasCleaning(Agent):
    """Syngas cooler (heat recovery), dust/tar removal (assumed complete), water quench to
    T_clean with condensation, and a blower that recovers the cleaning pressure drop."""

    name = "gas_cleaning"

    def __init__(self, params: CleaningParams = CleaningParams(), name: Optional[str] = None):
        super().__init__(name)
        self.p = params

    def process(self, gas: GasStream) -> Dict[str, Any]:
        p = self.p
        q_rec = p.heat_recovery_eff * sensible_heat(gas, gas.T, p.T_heat_recovery_min)
        low = gas.copy()
        low.p = gas.p - p.pressure_drop
        clean, cond = condense_water(low, p.T_clean)
        w_blower = compression_power(clean, gas.p, p.blower_eta)
        clean.p = gas.p
        return {"gas": clean, "heat_recovered_W": q_rec, "condensate_mol_s": cond,
                "electricity_W": w_blower}

    def step(self, t, dt, msg):
        out = self.process(msg["gas"])
        self.book(dt, electricity_W=out["electricity_W"], heat_recovered_W=out["heat_recovered_W"])
        return out


# ------------------------------------------------------------------------------------------
@dataclass
class MethanationParams:
    level: str = "none"  # "none" | "partial" | "full"
    h2_fraction_of_stoich: Optional[float] = None  # override: share of stoichiometric H2 make-up
    T: float = 300 + 273.15  # outlet of last (equilibrium) reactor stage
    p: float = 10e5
    compressor_eta: float = 0.75
    compressor_stages: int = 2

    @property
    def h2_share(self) -> float:
        if self.h2_fraction_of_stoich is not None:
            return self.h2_fraction_of_stoich
        return {"none": 0.0, "partial": 0.5, "full": 1.0}[self.level]


def methanation_equilibrium(n0: Dict[str, float], T: float, p: float) -> Dict[str, float]:
    """Chemical equilibrium of CO + 3H2 <-> CH4 + H2O and CO + H2O <-> CO2 + H2 at (T, p).
    Inerts (N2, Ar) are carried through. Solid carbon is not considered.
    Solved by nested bisection on the two reaction extents in log-residual form."""
    K1 = th.equilibrium_constant({"CO": -1, "H2": -3, "CH4": 1, "H2O": 1}, T)
    K2 = th.equilibrium_constant({"CO": -1, "H2O": -1, "CO2": 1, "H2": 1}, T)
    CO0, CO20, H20 = n0.get("CO", 0.0), n0.get("CO2", 0.0), n0.get("H2", 0.0)
    W0, M0 = n0.get("H2O", 0.0), n0.get("CH4", 0.0)
    inert = sum(v for k, v in n0.items() if k in ("N2", "Ar"))
    P = p / th.P_REF
    tiny = 1e-30

    def wgs(x1: float):
        lo = max(-CO20, 3 * x1 - H20)
        hi = min(CO0 - x1, W0 + x1)
        if lo > hi:
            return None
        def r(x2):
            co2, h2 = CO20 + x2, H20 - 3 * x1 + x2
            co, w = CO0 - x1 - x2, W0 + x1 - x2
            return (math.log(max(co2, tiny)) + math.log(max(h2, tiny))
                    - math.log(max(co, tiny)) - math.log(max(w, tiny)) - math.log(K2))
        a, b = lo, hi
        for _ in range(200):
            m = 0.5 * (a + b)
            if r(m) > 0:
                b = m
            else:
                a = m
        return 0.5 * (a + b)

    def comp(x1, x2):
        return {"CO": CO0 - x1 - x2, "CO2": CO20 + x2, "H2": H20 - 3 * x1 + x2,
                "H2O": W0 + x1 - x2, "CH4": M0 + x1}

    # upper bound of x1: largest value for which a feasible x2 exists
    x1_hi = (CO0 + CO20 + H20 + W0)  # loose; tighten by bisection on feasibility
    lo, hi = 0.0, x1_hi
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if wgs(m) is None:
            hi = m
        else:
            lo = m
    x1_max = lo

    def r1(x1):
        x2 = wgs(x1)
        n = comp(x1, x2)
        tot = sum(n.values()) + inert
        y = {k: max(v, tiny) / tot for k, v in n.items()}
        return (math.log(y["CH4"]) + math.log(y["H2O"]) - math.log(y["CO"])
                - 3 * math.log(y["H2"]) - 2 * math.log(P) - math.log(K1))

    a, b = (0.0, x1_max)
    if r1(a) >= 0:
        x1 = 0.0
    else:
        for _ in range(200):
            m = 0.5 * (a + b)
            if r1(m) > 0:
                b = m
            else:
                a = m
        x1 = 0.5 * (a + b)
    out = comp(x1, wgs(x1))
    out = {k: max(v, 0.0) for k, v in out.items()}
    for k in ("N2", "Ar"):
        if n0.get(k):
            out[k] = n0[k]
    return out


class Methanation(Agent):
    """Compression to methanation pressure, H2 make-up, equilibrium methanation at (T, p),
    reaction-heat recovery and final drying to the cleaning temperature."""

    name = "methanation"

    def __init__(self, params: MethanationParams = MethanationParams(),
                 T_dry: float = 40 + 273.15, name: Optional[str] = None):
        super().__init__(name)
        self.p = params
        self.T_dry = T_dry

    def h2_demand(self, gas: GasStream) -> float:
        """External H2 (mol/s) for the configured share of the stoichiometric make-up to
        module M = (H2 - CO2)/(CO + CO2) = 3."""
        n = gas.n
        need = 3 * (n.get("CO", 0) + n.get("CO2", 0)) + n.get("CO2", 0) - n.get("H2", 0)
        return max(0.0, need) * self.p.h2_share

    def process(self, gas: GasStream, h2_mol_s: Optional[float] = None) -> Dict[str, Any]:
        p = self.p
        if p.level == "none":
            return {"gas": gas, "h2_mol_s": 0.0, "electricity_W": 0.0, "heat_recovered_W": 0.0}
        h2 = self.h2_demand(gas) if h2_mol_s is None else h2_mol_s
        w_comp = compression_power(gas, p.p, p.compressor_eta, p.compressor_stages)
        feed = dict(gas.n)
        feed["H2"] = feed.get("H2", 0.0) + h2
        H_in = th.mix_h(feed, gas.T)  # H2 assumed delivered at the same temperature
        eq = methanation_equilibrium(feed, p.T, p.p)
        q_reac = H_in - th.mix_h(eq, p.T)  # heat removed from reactors (exothermic > 0)
        out = GasStream(eq, p.T, p.p)
        q_cool = sensible_heat(out, p.T, self.T_dry)
        dried, cond = condense_water(out, self.T_dry)
        return {"gas": dried, "h2_mol_s": h2, "electricity_W": w_comp,
                "heat_recovered_W": q_reac, "low_grade_heat_W": q_cool, "condensate_mol_s": cond}

    def step(self, t, dt, msg):
        out = self.process(msg["gas"], msg.get("h2_mol_s"))
        self.book(dt, h2_mol_s=out["h2_mol_s"], electricity_W=out["electricity_W"])
        return out


# ------------------------------------------------------------------------------------------
@dataclass
class ElectrolyserParams:
    specific_energy_kwh_per_kg: float = 55.0  # system, incl. BoP (assumption, see docs)
    capacity_W: float = math.inf  # electric capacity; inf = sized to demand
    min_load: float = 0.0  # fraction of capacity
    water_kg_per_kg_h2: float = 9.0  # stoichiometric 8.94 kg/kg


class Electrolyser(Agent):
    """Water electrolysis: electricity -> H2 + 1/2 O2 (stoichiometric O2 co-product)."""

    name = "electrolyser"

    def __init__(self, params: ElectrolyserParams = ElectrolyserParams(),
                 name: Optional[str] = None):
        super().__init__(name)
        self.p = params

    def power_for_h2(self, h2_mol_s: float) -> float:
        kg_s = h2_mol_s * th.MOLAR_MASS["H2"] / 1000
        return kg_s * self.p.specific_energy_kwh_per_kg * KWH

    def process(self, power_W: float) -> Dict[str, float]:
        cap = self.p.capacity_W
        if power_W < 0:
            raise ValueError("power must be >= 0")
        if power_W > 0 and math.isfinite(cap):
            power_W = min(power_W, cap)
            if power_W < self.p.min_load * cap:
                power_W = 0.0
        h2_kg_s = power_W / (self.p.specific_energy_kwh_per_kg * KWH)
        h2 = h2_kg_s * 1000 / th.MOLAR_MASS["H2"]
        return {"electricity_W": power_W, "h2_mol_s": h2, "o2_mol_s": h2 / 2,
                "water_kg_s": h2_kg_s * self.p.water_kg_per_kg_h2}

    def step(self, t, dt, msg):
        out = self.process(msg["power_W"])
        self.book(dt, **{k: out[k] for k in ("electricity_W", "h2_mol_s", "o2_mol_s")})
        return out


# ------------------------------------------------------------------------------------------
class GasStorage(Agent):
    """Ideal gas buffer (H2 or O2) with capacity in mol, charge/discharge rate limits and a
    specific compression energy for charging (W per mol/s)."""

    def __init__(self, name: str, capacity_mol: float = 0.0, max_rate_mol_s: float = math.inf,
                 soc0: float = 0.0, charge_energy_J_per_mol: float = 0.0):
        super().__init__(name)
        self.capacity = capacity_mol
        self.max_rate = max_rate_mol_s
        self.level = soc0 * capacity_mol
        self.level0 = self.level
        self.e_charge = charge_energy_J_per_mol

    def reset(self) -> None:
        super().reset()
        self.level = self.level0

    def exchange(self, net_in_mol_s: float, dt: float) -> Dict[str, float]:
        """Try to absorb (>0) or deliver (<0) ``net_in_mol_s`` for dt seconds. Returns the
        realised flow and the part that could not be handled (spill > 0 / shortfall < 0)."""
        rate = max(-self.max_rate, min(self.max_rate, net_in_mol_s))
        if rate > 0:
            rate = min(rate, (self.capacity - self.level) / dt)
        else:
            rate = max(rate, -self.level / dt)
        self.level += rate * dt
        return {"flow_mol_s": rate, "unserved_mol_s": net_in_mol_s - rate,
                "electricity_W": max(rate, 0.0) * self.e_charge, "level_mol": self.level}

    def step(self, t, dt, msg):
        out = self.exchange(msg["net_in_mol_s"], dt)
        self.book(dt, electricity_W=out["electricity_W"], flow_mol_s=out["flow_mol_s"])
        self.last["level_mol"] = self.level
        return out


# ------------------------------------------------------------------------------------------
@dataclass
class ASUParams:
    specific_energy_kwh_per_kg_o2: float = 0.25  # cryogenic ASU, assumption (see docs)


class AirSeparationUnit(Agent):
    """Oxygen supply from a cryogenic air separation unit (electricity only)."""

    name = "asu"

    def __init__(self, params: ASUParams = ASUParams(), name: Optional[str] = None):
        super().__init__(name)
        self.p = params

    def process(self, o2_mol_s: float) -> Dict[str, float]:
        kg_s = o2_mol_s * th.MOLAR_MASS["O2"] / 1000
        return {"o2_mol_s": o2_mol_s, "electricity_W": kg_s * self.p.specific_energy_kwh_per_kg_o2 * KWH}

    def step(self, t, dt, msg):
        out = self.process(msg["o2_mol_s"])
        self.book(dt, **out)
        return out


# ------------------------------------------------------------------------------------------
@dataclass
class BurnerParams:
    excess_air: float = 1.10  # lambda
    T_air: float = th.T_REF  # combustion-air temperature (recuperative preheat if > 25 degC)
    T_fuel: float = 40 + 273.15


class Burner(Agent):
    """Industrial burner / furnace. Complete combustion with air at lambda; the useful
    (available) heat is the enthalpy released by cooling the flue gas from the adiabatic
    flame temperature down to the flue-gas exit temperature of the process."""

    name = "burner"

    def __init__(self, params: BurnerParams = BurnerParams(), name: Optional[str] = None):
        super().__init__(name)
        self.p = params

    def flue_gas(self, fuel: GasStream) -> tuple[Dict[str, float], Dict[str, float]]:
        at = fuel.atoms()
        o2_stoich = at["C"] + at["H"] / 4 - at["O"] / 2
        o2 = max(o2_stoich, 0.0) * self.p.excess_air
        air = {k: o2 * v / th.AIR["O2"] for k, v in th.AIR.items()}
        flue = {"CO2": at["C"], "H2O": at["H"] / 2, "N2": at["N"] / 2 + air["N2"],
                "O2": o2 - o2_stoich, "Ar": air["Ar"] + fuel.n.get("Ar", 0.0)}
        return air, flue

    def adiabatic_flame_temperature(self, fuel: GasStream) -> float:
        air, flue = self.flue_gas(fuel)
        H = th.mix_h(fuel.n, self.p.T_fuel) + th.mix_h(air, self.p.T_air)
        return th.solve_temperature(flue, H)

    def available_heat(self, fuel: GasStream, T_flue_exit: float) -> float:
        """Useful heat (W) delivered above T_flue_exit."""
        air, flue = self.flue_gas(fuel)
        H = th.mix_h(fuel.n, self.p.T_fuel) + th.mix_h(air, self.p.T_air)
        return H - th.mix_h(flue, T_flue_exit)

    def available_heat_fraction(self, fuel: GasStream, T_flue_exit: float) -> float:
        return self.available_heat(fuel, T_flue_exit) / fuel.lhv_power

    def co2_kg_s(self, fuel: GasStream) -> float:
        """CO2 in the flue gas (kg/s); fossil for natural gas, biogenic for syngas."""
        return fuel.atoms()["C"] * th.MOLAR_MASS["CO2"] / 1000

    def step(self, t, dt, msg):
        q = self.available_heat(msg["fuel"], msg["T_flue_exit"])
        self.book(dt, heat_W=q, fuel_lhv_W=msg["fuel"].lhv_power)
        return {"heat_W": q}


# ------------------------------------------------------------------------------------------
class ProcessHeatDemand(Agent):
    """Industrial process requesting high-temperature heat (W) at a process temperature."""

    name = "process_heat_demand"

    def __init__(self, profile_W, T_process: float = 1000 + 273.15,
                 flue_exit_approach_K: float = 0.0, name: Optional[str] = None):
        super().__init__(name)
        self.profile = list(profile_W)
        self.T_process = T_process
        self.approach = flue_exit_approach_K

    @property
    def T_flue_exit(self) -> float:
        return self.T_process + self.approach

    def step(self, t, dt, msg):
        q = self.profile[t % len(self.profile)]
        self.book(dt, demand_W=q)
        return {"demand_W": q, "T_flue_exit": self.T_flue_exit}


class GridMarket(Agent):
    """Electricity market: hourly price (EUR/MWh) and the plant's cumulative purchases."""

    name = "grid_market"

    def __init__(self, prices_eur_mwh=None, name: Optional[str] = None):
        super().__init__(name)
        self.prices = list(prices_eur_mwh) if prices_eur_mwh is not None else [0.0]

    def price(self, t: int) -> float:
        return self.prices[t % len(self.prices)]

    def step(self, t, dt, msg):
        p = self.price(t)
        w = msg["electricity_W"]
        self.book(dt, electricity_W=w, cost_eur_s=w / MW * p / 3600)
        return {"price": p}
