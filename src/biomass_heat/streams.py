"""Gas stream container used to pass material between agents."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

from . import thermo as th


@dataclass
class GasStream:
    """Molar flows (mol/s) of ideal-gas species at temperature T (K) and pressure p (Pa)."""

    n: Dict[str, float] = field(default_factory=dict)
    T: float = th.T_REF
    p: float = th.P_REF

    def copy(self) -> "GasStream":
        return GasStream(dict(self.n), self.T, self.p)

    # --- totals ---------------------------------------------------------------------------
    @property
    def total(self) -> float:
        return sum(self.n.values())

    @property
    def dry_total(self) -> float:
        return self.total - self.n.get("H2O", 0.0)

    def x(self, dry: bool = False) -> Dict[str, float]:
        """Mole fractions (optionally on a dry basis)."""
        tot = self.dry_total if dry else self.total
        return {k: v / tot for k, v in self.n.items() if v > 0 and not (dry and k == "H2O")}

    @property
    def mass_flow(self) -> float:
        """kg/s"""
        return sum(v * th.MOLAR_MASS[k] for k, v in self.n.items()) / 1000.0

    @property
    def molar_mass(self) -> float:
        return 1000.0 * self.mass_flow / self.total

    # --- energy ---------------------------------------------------------------------------
    @property
    def enthalpy(self) -> float:
        """Absolute enthalpy flow incl. formation, W."""
        return th.mix_h(self.n, self.T)

    @property
    def lhv_power(self) -> float:
        """Chemical energy flow on LHV basis, W."""
        return sum(v * th.lhv_molar(k) for k, v in self.n.items())

    @property
    def hhv_power(self) -> float:
        return sum(v * th.hhv_molar(k) for k, v in self.n.items())

    def atoms(self) -> Dict[str, float]:
        out = {"C": 0.0, "H": 0.0, "O": 0.0, "N": 0.0}
        for k, v in self.n.items():
            c, h, o, nn = th.ATOMS[k]
            out["C"] += c * v
            out["H"] += h * v
            out["O"] += o * v
            out["N"] += nn * v
        return out

    # --- fuel-gas quality (dry basis, normal conditions 0 degC / 1.01325 bar) ----------------
    def lhv_volumetric(self, dry: bool = True) -> float:
        """MJ/Nm3"""
        tot = self.dry_total if dry else self.total
        return self.lhv_power / (tot * th.V_MOLAR_NORMAL) / 1e6

    def hhv_volumetric(self, dry: bool = True) -> float:
        """MJ/Nm3"""
        tot = self.dry_total if dry else self.total
        return self.hhv_power / (tot * th.V_MOLAR_NORMAL) / 1e6

    def relative_density(self, dry: bool = True) -> float:
        n = {k: v for k, v in self.n.items() if not (dry and k == "H2O")}
        m = sum(v * th.MOLAR_MASS[k] for k, v in n.items()) / sum(n.values())
        return m / th.M_AIR

    def wobbe_index(self, dry: bool = True) -> float:
        """Upper (HHV-based) Wobbe index, MJ/Nm3."""
        return self.hhv_volumetric(dry) / self.relative_density(dry) ** 0.5

    def __add__(self, other: "GasStream") -> "GasStream":
        """Adiabatic mixing (enthalpy conserving); pressure = min of both."""
        n = dict(self.n)
        for k, v in other.n.items():
            n[k] = n.get(k, 0.0) + v
        H = self.enthalpy + other.enthalpy
        out = GasStream(n, self.T, min(self.p, other.p))
        if out.total > 0:
            out.T = th.solve_temperature(n, H)
        return out

    def scaled(self, f: float) -> "GasStream":
        return GasStream({k: v * f for k, v in self.n.items()}, self.T, self.p)
