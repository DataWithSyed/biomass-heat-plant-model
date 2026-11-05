"""Run every plant configuration at the reference design point and write results/ tables
and figures. Usage:  python scripts/run_configurations.py"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from biomass_heat import ALL_CONFIGS, CONFIGS, Burner, Plant, natural_gas

OUT = Path(__file__).resolve().parents[1] / "results"
HEAT_MW = 10.0
T_PROCESS_C = 1000.0
SHORT = {"air": "Air", "oxy_none": "Oxy (ASU)", "oxy_partial": "Oxy + partial meth.",
         "oxy_full": "Oxy + full meth.", "oxy_none_elec": "Oxy (electrolyser O2, H2 blend)",
         "steam_dfb": "Steam (DFB)"}


def to_markdown(df: pd.DataFrame) -> str:
    """Minimal Markdown table writer (avoids an extra dependency)."""
    fmt = lambda v: "" if pd.isna(v) else (
        (f"{v:.0f}" if abs(v) >= 100 else f"{v:.3g}") if isinstance(v, float) else str(v))
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines) + "\n"


def main() -> pd.DataFrame:
    OUT.mkdir(exist_ok=True)
    rows = [Plant(k).design_point(HEAT_MW * 1e6, T_PROCESS_C + 273.15) for k in ALL_CONFIGS]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "design_point_all_kpis.csv", index=False, float_format="%.5g")

    ng = natural_gas()
    burner = Burner()
    ng_row = {"config": "natural_gas", "fuel_lhv_MJ_Nm3": ng.lhv_volumetric(),
              "fuel_wobbe_MJ_Nm3": ng.wobbe_index(),
              "adiabatic_flame_C": burner.adiabatic_flame_temperature(ng) - 273.15,
              "available_heat_frac": burner.available_heat_fraction(ng, T_PROCESS_C + 273.15)}

    cols = {"config": "Configuration", "cold_gas_efficiency": "CGE (-)",
            "fuel_lhv_MJ_Nm3": "Fuel LHV (MJ/Nm3)", "fuel_wobbe_MJ_Nm3": "Wobbe (MJ/Nm3)",
            "adiabatic_flame_C": "T_ad (degC)", "available_heat_frac": "Avail. heat @1000 degC",
            "biomass_MWh_per_MWh": "Biomass MWh/MWh_heat",
            "electricity_MWh_per_MWh": "Electricity MWh/MWh_heat",
            "o2_kg_per_MWh": "O2 kg/MWh_heat", "h2_kg_per_MWh": "H2 kg/MWh_heat",
            "o2_surplus_kg_h": "O2 surplus (kg/h)", "chain_eff_heat": "Heat / (biomass+el.)"}
    tab = df[list(cols)].copy()
    tab = pd.concat([tab, pd.DataFrame([ng_row])], ignore_index=True)
    tab = tab.rename(columns=cols)
    tab.to_csv(OUT / "results_table.csv", index=False, float_format="%.4g")
    with open(OUT / "results_table.md", "w") as f:
        f.write(to_markdown(tab))

    comp = df.set_index("config")[[c for c in df.columns if c.startswith("fuel_") and "vol" in c]]
    comp.round(2).to_csv(OUT / "fuel_gas_composition.csv")

    # ---- figure 1: fuel quality vs natural gas ---------------------------------------
    keys = list(ALL_CONFIGS)
    lab = [SHORT[k] for k in keys]
    d = df.set_index("config").loc[keys]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    x = np.arange(len(keys))
    for a, col, title, ref in [
        (ax[0], "fuel_lhv_MJ_Nm3", "Fuel-gas LHV (dry, MJ/Nm$^3$)", ng.lhv_volumetric()),
        (ax[1], "fuel_wobbe_MJ_Nm3", "Upper Wobbe index (MJ/Nm$^3$)", ng.wobbe_index()),
        (ax[2], "adiabatic_flame_C", "Adiabatic flame temp. (°C, λ=1.1, cold air)",
         ng_row["adiabatic_flame_C"])]:
        a.bar(x, d[col], color=["#8c6d31" if k in CONFIGS else "#bdbdbd" for k in keys])
        a.axhline(ref, color="k", ls="--", lw=1)
        a.text(-0.4, ref, "natural gas", va="bottom", ha="left", fontsize=8)
        a.set_ylim(0, max(ref, d[col].max()) * 1.12)
        a.set_xticks(x, lab, rotation=35, ha="right", fontsize=8)
        a.set_title(title, fontsize=10)
    fig.suptitle("Fuel-gas quality of biomass process chains vs natural gas "
                 "(grey = extra variants)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fuel_quality_vs_natural_gas.png", dpi=150)
    plt.close(fig)

    # ---- figure 2: inputs per MWh of process heat --------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].bar(x, d["biomass_MWh_per_MWh"], label="biomass (LHV)", color="#8c6d31")
    ax[0].bar(x, d["electricity_MWh_per_MWh"], bottom=d["biomass_MWh_per_MWh"],
              label="electricity", color="#3182bd")
    ax[0].axhline(d["ng_MWh_per_MWh"].iloc[0], color="k", ls="--", lw=1,
                  label="natural gas (LHV) reference")
    ax[0].set_ylabel("MWh input per MWh process heat @1000 °C")
    ax[0].legend(fontsize=8)
    w = 0.38
    ax[1].bar(x - w / 2, d["o2_kg_per_MWh"], w, label="O$_2$ demand of gasifier", color="#e6550d")
    ax[1].bar(x + w / 2, d["o2_coproduct_kg_h"] / HEAT_MW, w, label="O$_2$ co-produced by electrolyser",
              color="#fdae6b")
    ax[1].plot(x, d["h2_kg_per_MWh"] * 10, "ks", label="H$_2$ demand ×10")
    ax[1].set_ylabel("kg per MWh process heat")
    ax[1].legend(fontsize=8)
    for a in ax:
        a.set_xticks(x, lab, rotation=35, ha="right", fontsize=8)
    fig.suptitle("Resource demand per MWh of high-temperature process heat", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "inputs_per_mwh_heat.png", dpi=150)
    plt.close(fig)

    # ---- figure 3: available heat vs process temperature -------------------------------
    temps = np.arange(200, 1501, 50)
    fig, a = plt.subplots(figsize=(7.5, 4.8))
    ah = {}
    for k in keys:
        fuel = Plant(k).chain(1.0)["fuel"]
        ah[k] = [burner.available_heat_fraction(fuel, t + 273.15) for t in temps]
        a.plot(temps, ah[k], label=SHORT[k], lw=2 if k in CONFIGS else 1,
               ls="-" if k in CONFIGS else ":")
    ah["natural_gas"] = [burner.available_heat_fraction(ng, t + 273.15) for t in temps]
    a.plot(temps, ah["natural_gas"], "k--", lw=2, label="natural gas")
    a.set_xlabel("Flue-gas exit temperature = process temperature (°C)")
    a.set_ylabel("Available heat / fuel LHV (-)")
    a.set_ylim(0, 1.05)
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
    a.set_title("Available heat of the fuel gases (λ=1.1, no air preheat)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "available_heat_vs_process_temperature.png", dpi=150)
    plt.close(fig)
    pd.DataFrame(ah, index=temps).rename_axis("T_process_C").round(4).to_csv(
        OUT / "available_heat_vs_temperature.csv")
    return tab


if __name__ == "__main__":
    print(main().to_string(index=False))
