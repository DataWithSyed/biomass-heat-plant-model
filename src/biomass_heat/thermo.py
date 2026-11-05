"""Ideal-gas thermodynamics from NASA 7-coefficient polynomials.

Coefficients are taken from the GRI-Mech 3.0 thermodynamic data set
(Smith et al., http://combustion.berkeley.edu/gri-mech/), as distributed with
Cantera's ``gri30.yaml``. Valid ~300-3500 K. All enthalpies are absolute
(include the enthalpy of formation at 298.15 K), in J/mol.
"""
from __future__ import annotations

import math
from typing import Mapping

R = 8.314462618  # J/(mol K)
T_REF = 298.15  # K
P_REF = 101325.0  # Pa
V_MOLAR_NORMAL = 0.022414  # m3/mol at 0 degC, 1.01325 bar (ideal gas, "Nm3")
L_H2O_25C = 44.0e3  # J/mol, latent heat of water at 25 degC (approx. 44.0 kJ/mol)

# species: (T_mid, low-T coeffs a1..a7, high-T coeffs a1..a7)
NASA7: dict[str, tuple[float, tuple[float, ...], tuple[float, ...]]] = {
    "H2": (1000.0, (2.34433112, 0.00798052075, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238), (3.3372792, -4.94024731e-05, 4.99456778e-07, -1.79566394e-10, 2.00255376e-14, -950.158922, -3.20502331)),
    "CO": (1000.0, (3.57953347, -0.00061035368, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -14344.086, 3.50840928), (2.71518561, 0.00206252743, -9.98825771e-07, 2.30053008e-10, -2.03647716e-14, -14151.8724, 7.81868772)),
    "CO2": (1000.0, (2.35677352, 0.00898459677, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -48371.9697, 9.90105222), (3.85746029, 0.00441437026, -2.21481404e-06, 5.23490188e-10, -4.72084164e-14, -48759.166, 2.27163806)),
    "CH4": (1000.0, (5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11, -10246.6476, -4.64130376), (0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -1.0181523e-13, -9468.34459, 18.437318)),
    "H2O": (1000.0, (4.19864056, -0.0020364341, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -30293.7267, -0.849032208), (3.03399249, 0.00217691804, -1.64072518e-07, -9.7041987e-11, 1.68200992e-14, -30004.2971, 4.9667701)),
    "N2": (1000.0, (3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372), (2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528)),
    "O2": (1000.0, (3.78245636, -0.00299673416, 9.84730201e-06, -9.68129509e-09, 3.24372837e-12, -1063.94356, 3.65767573), (3.28253784, 0.00148308754, -7.57966669e-07, 2.09470555e-10, -2.16717794e-14, -1088.45772, 5.45323129)),
    "Ar": (1000.0, (2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366), (2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366)),
    "C2H6": (1000.0, (4.29142492, -0.0055015427, 5.99438288e-05, -7.08466285e-08, 2.68685771e-11, -11522.2055, 2.66682316), (1.0718815, 0.0216852677, -1.00256067e-05, 2.21412001e-09, -1.9000289e-13, -11426.3932, 15.1156107)),
    "C3H8": (1000.0, (0.93355381, 0.026424579, 6.1059727e-06, -2.1977499e-08, 9.5149253e-12, -13958.52, 19.201691), (7.5341368, 0.018872239, -6.2718491e-06, 9.1475649e-10, -4.7838069e-14, -16467.516, -17.892349)),}

MOLAR_MASS = {  # g/mol
    "H2": 2.016, "CO": 28.010, "CO2": 44.009, "CH4": 16.043, "H2O": 18.015,
    "N2": 28.014, "O2": 31.998, "Ar": 39.948, "C2H6": 30.069, "C3H8": 44.096,
}
ATOMS = {  # (C, H, O, N)
    "H2": (0, 2, 0, 0), "CO": (1, 0, 1, 0), "CO2": (1, 0, 2, 0), "CH4": (1, 4, 0, 0),
    "H2O": (0, 2, 1, 0), "N2": (0, 0, 0, 2), "O2": (0, 0, 2, 0), "Ar": (0, 0, 0, 0),
    "C2H6": (2, 6, 0, 0), "C3H8": (3, 8, 0, 0),
}
SPECIES = tuple(NASA7)
M_AIR = 28.965  # g/mol, dry air
AIR = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0097}  # dry air mole fractions (sum 1.0)


def _coeffs(sp: str, T: float) -> tuple[float, ...]:
    tmid, low, high = NASA7[sp]
    return low if T < tmid else high


def cp(sp: str, T: float) -> float:
    """Molar heat capacity, J/(mol K)."""
    a = _coeffs(sp, T)
    return R * (a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3 + a[4] * T**4)


def h(sp: str, T: float) -> float:
    """Absolute molar enthalpy (incl. formation), J/mol."""
    a = _coeffs(sp, T)
    return R * T * (a[0] + a[1] * T / 2 + a[2] * T**2 / 3 + a[3] * T**3 / 4
                    + a[4] * T**4 / 5 + a[5] / T)


def s0(sp: str, T: float) -> float:
    """Standard-state molar entropy at 1 atm, J/(mol K)."""
    a = _coeffs(sp, T)
    return R * (a[0] * math.log(T) + a[1] * T + a[2] * T**2 / 2 + a[3] * T**3 / 3
                + a[4] * T**4 / 4 + a[6])


def g0(sp: str, T: float) -> float:
    """Standard-state molar Gibbs energy, J/mol."""
    return h(sp, T) - T * s0(sp, T)


def hf(sp: str) -> float:
    """Enthalpy of formation at 298.15 K, J/mol."""
    return h(sp, T_REF)


def equilibrium_constant(nu: Mapping[str, float], T: float) -> float:
    """Kp (pressure in atm/bar-ish standard state) for reaction with stoich. coefficients nu
    (products positive, reactants negative)."""
    dg = sum(v * g0(sp, T) for sp, v in nu.items())
    return math.exp(-dg / (R * T))


def lhv_molar(sp: str) -> float:
    """Lower heating value at 25 degC, J/mol (0 for inert species)."""
    c, hh, o, _ = ATOMS[sp]
    if sp in ("CO2", "H2O", "N2", "O2", "Ar"):
        return 0.0
    o2 = c + hh / 4 - o / 2
    return h(sp, T_REF) + o2 * h("O2", T_REF) - c * h("CO2", T_REF) - hh / 2 * h("H2O", T_REF)


def hhv_molar(sp: str) -> float:
    """Higher heating value at 25 degC, J/mol (water condensed)."""
    hh = ATOMS[sp][1]
    return lhv_molar(sp) + (hh / 2 * L_H2O_25C if lhv_molar(sp) > 0 else 0.0)


def mix_h(n: Mapping[str, float], T: float) -> float:
    """Total enthalpy of a gas mixture with molar amounts/flows n at T, J (or W if n in mol/s)."""
    return sum(v * h(sp, T) for sp, v in n.items() if v)


def solve_temperature(n: Mapping[str, float], H_target: float, lo: float = 250.0,
                      hi: float = 3500.0) -> float:
    """Temperature at which mixture n has absolute enthalpy H_target (bisection)."""
    f_lo, f_hi = mix_h(n, lo) - H_target, mix_h(n, hi) - H_target
    if f_lo * f_hi > 0:
        raise ValueError("target enthalpy outside temperature bracket")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f = mix_h(n, mid) - H_target
        if abs(f) < 1e-9 * max(1.0, abs(H_target)) or hi - lo < 1e-6:
            return mid
        if (f < 0) == (f_lo < 0):
            lo, f_lo = mid, f
        else:
            hi = mid
    return 0.5 * (lo + hi)


def p_sat_water(T: float) -> float:
    """Saturation pressure of water over liquid, Pa (Arden Buck equation, 1996 revision of Buck 1981).
    """
    tc = T - 273.15
    return 611.21 * math.exp((18.678 - tc / 234.5) * (tc / (257.14 + tc)))
