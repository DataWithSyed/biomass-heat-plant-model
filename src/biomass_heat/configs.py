"""Pre-defined process-chain configurations.

All numbers are *engineering assumptions* for a screening study (see docs/assumptions.md);
none of them is a measured value from a specific plant.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict

from .agents import (ASUParams, BurnerParams, CleaningParams, ElectrolyserParams,
                     GasifierParams, MethanationParams)
from .biomass import Biomass
from .streams import GasStream


@dataclass(frozen=True)
class PlantConfig:
    key: str
    label: str
    gasifier: GasifierParams
    methanation: MethanationParams = field(default_factory=MethanationParams)
    o2_supply: str = "none"  # "none" | "asu" | "electrolyser"
    blend_surplus_h2: bool = False  # electrolyser H2 not used by methanation goes to fuel gas
    biomass: Biomass = field(default_factory=Biomass)
    cleaning: CleaningParams = field(default_factory=CleaningParams)
    electrolyser: ElectrolyserParams = field(default_factory=ElectrolyserParams)
    asu: ASUParams = field(default_factory=ASUParams)
    burner: BurnerParams = field(default_factory=BurnerParams)

    def with_(self, **kw) -> "PlantConfig":
        return replace(self, **kw)


AIR = GasifierParams(mode="air", carbon_conversion=0.95, ch4_carbon_fraction=0.08)
OXY = GasifierParams(mode="oxyfuel", carbon_conversion=0.97, ch4_carbon_fraction=0.10,
                     steam_to_dry_biomass=0.5)
DFB = GasifierParams(mode="steam", carbon_conversion=0.75, ch4_carbon_fraction=0.12,
                     steam_to_dry_biomass=0.75)

CONFIGS: Dict[str, PlantConfig] = {
    "air": PlantConfig("air", "Air gasification", AIR),
    "oxy_none": PlantConfig("oxy_none", "Oxyfuel, no methanation (O2 from ASU)", OXY,
                            o2_supply="asu"),
    "oxy_partial": PlantConfig("oxy_partial", "Oxyfuel + partial methanation (50 % H2)", OXY,
                               MethanationParams(level="partial"), o2_supply="electrolyser"),
    "oxy_full": PlantConfig("oxy_full", "Oxyfuel + full methanation", OXY,
                            MethanationParams(level="full"), o2_supply="electrolyser"),
}

EXTRA_CONFIGS: Dict[str, PlantConfig] = {
    "oxy_none_elec": PlantConfig("oxy_none_elec",
                                 "Oxyfuel, no methanation (O2 from electrolyser, H2 blended)",
                                 OXY, o2_supply="electrolyser", blend_surplus_h2=True),
    "steam_dfb": PlantConfig("steam_dfb", "Allothermal steam gasification (DFB)", DFB),
}

ALL_CONFIGS: Dict[str, PlantConfig] = {**CONFIGS, **EXTRA_CONFIGS}

# Illustrative natural gas (H-gas) composition, mol-% - an assumption, not a measured analysis.
NATURAL_GAS_X = {"CH4": 0.93, "C2H6": 0.04, "C3H8": 0.01, "CO2": 0.01, "N2": 0.01}


def natural_gas(mol_s: float = 1.0) -> GasStream:
    return GasStream({k: v * mol_s for k, v in NATURAL_GAS_X.items()}, 298.15)
