"""Small generic unit-operation helpers shared by several agents."""
from __future__ import annotations

from . import thermo as th
from .streams import GasStream


def compression_power(stream: GasStream, p_out: float, eta_is: float = 0.75,
                      stages: int = 1, T_intercool: float | None = None) -> float:
    """Electric/shaft power (W) to compress an ideal-gas stream from stream.p to p_out with
    ``stages`` equal-pressure-ratio stages and intercooling back to the inlet temperature
    (or ``T_intercool``). Uses cp at stage inlet temperature (cold-air-standard style)."""
    if p_out <= stream.p or stream.total <= 0:
        return 0.0
    T1 = T_intercool if T_intercool is not None else stream.T
    pr = (p_out / stream.p) ** (1.0 / stages)
    n_tot = stream.total
    cp_mix = sum(v * th.cp(k, T1) for k, v in stream.n.items()) / n_tot
    exponent = th.R / cp_mix
    w_stage = n_tot * cp_mix * T1 * (pr ** exponent - 1.0) / eta_is
    return stages * w_stage


def condense_water(stream: GasStream, T_out: float) -> tuple[GasStream, float]:
    """Cool a stream to T_out and condense water down to saturation at stream.p.
    Returns (gas, condensed water mol/s)."""
    out = stream.copy()
    out.T = T_out
    w = out.n.get("H2O", 0.0)
    dry = out.dry_total
    y_sat = min(th.p_sat_water(T_out) / out.p, 0.999)
    w_max = y_sat * dry / (1 - y_sat)
    cond = max(0.0, w - w_max)
    out.n["H2O"] = w - cond
    return out, cond


def sensible_heat(stream: GasStream, T_from: float, T_to: float) -> float:
    """Enthalpy released (W) when cooling ``stream`` from T_from to T_to, no condensation."""
    return th.mix_h(stream.n, T_from) - th.mix_h(stream.n, T_to)
