"""Solid biomass fuel description (ultimate analysis) and heating values."""
from __future__ import annotations

from dataclasses import dataclass

from . import thermo as th

M = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "S": 32.06}


@dataclass(frozen=True)
class Biomass:
    """Ultimate analysis in wt-% on a dry basis (must sum to 100 incl. ash).

    The default is an *assumed* generic woody biomass (wood chips) whose composition lies
    within the range of woody-biomass entries of the Phyllis2 database (TNO/ECN,
    https://phyllis.nl). It is not a measured sample.
    """

    name: str = "wood chips (assumed generic woody biomass)"
    C: float = 50.0
    H: float = 6.0
    O: float = 43.3
    N: float = 0.2
    S: float = 0.0
    ash: float = 0.5
    moisture: float = 15.0  # wt-% as received (wet basis), assumption: dried chips

    def __post_init__(self) -> None:
        tot = self.C + self.H + self.O + self.N + self.S + self.ash
        if abs(tot - 100.0) > 0.05:
            raise ValueError(f"ultimate analysis must sum to 100 wt-% (got {tot:.2f})")
        if not 0 <= self.moisture < 100:
            raise ValueError("moisture must be in [0, 100)")

    # --- heating values ---------------------------------------------------------------
    @property
    def hhv_dry(self) -> float:
        """HHV of dry biomass, MJ/kg, Channiwala & Parikh (2002) unified correlation."""
        return (0.3491 * self.C + 1.1783 * self.H + 0.1005 * self.S - 0.1034 * self.O
                - 0.0151 * self.N - 0.0211 * self.ash)

    @property
    def lhv_dry(self) -> float:
        """LHV of dry biomass, MJ/kg (subtract latent heat of the water formed from H)."""
        water_kg = self.H / 100 * 18.015 / (2 * 1.008)
        return self.hhv_dry - water_kg * th.L_H2O_25C / 18.015 / 1000

    @property
    def lhv_ar(self) -> float:
        """LHV as received, MJ/kg (moisture must be evaporated)."""
        w = self.moisture / 100
        return self.lhv_dry * (1 - w) - w * th.L_H2O_25C / 18.015 / 1000

    def atoms_per_kg_dry(self) -> dict[str, float]:
        """mol of element per kg of dry biomass."""
        return {e: getattr(self, e) * 10.0 / M[e] for e in ("C", "H", "O", "N")}

    def formation_enthalpy_dry(self) -> float:
        """Enthalpy of formation of 1 kg dry biomass at 25 degC, J/kg, from the HHV
        (complete combustion to CO2(g), liquid H2O and N2)."""
        a = self.atoms_per_kg_dry()
        h_h2o_liq = th.hf("H2O") - th.L_H2O_25C
        return self.hhv_dry * 1e6 + a["C"] * th.hf("CO2") + a["H"] / 2 * h_h2o_liq
