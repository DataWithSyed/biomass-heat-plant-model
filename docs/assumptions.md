# Assumption register

Every number below is an **engineering assumption for a screening model**. None is a
measured value from a specific plant. "Source" says where the value or its range comes from;
"assumption" means a round value chosen by the author within ranges reported in the cited
reviews and to be varied in the sensitivity study (`biomass-heat-lcoh`).

| # | Item | Value | Label | Source / rationale |
|---|---|---|---|---|
| A1 | Biomass ultimate analysis (dry, wt-%) | C 50.0, H 6.0, O 43.3, N 0.2, ash 0.5 | assumption | Generic woody biomass, within woody entries of the Phyllis2 database (TNO, phyllis.nl). Not a measured sample. |
| A2 | Biomass moisture (as received) | 15 wt-% | assumption | Dried wood chips. |
| A3 | Biomass HHV | Channiwala & Parikh (2002) correlation -> 20.0 MJ/kg dry | correlation | Channiwala & Parikh, Fuel 81 (2002) 1051-1063, doi:10.1016/S0016-2361(01)00131-4 |
| A4 | Gasifier outlet temperature | 850 degC (all modes) | assumption | Typical fluidised-bed range 800-900 degC (Molino et al. 2016; Basu 2018). |
| A5 | Carbon conversion | air 0.95, oxyfuel 0.97, DFB gasification reactor 0.75 (char to combustor) | assumption | Screening values; DFB char split is what heats the combustor. |
| A6 | CH4 yield (fraction of biomass C) | air 0.08, oxyfuel 0.10, DFB 0.12 | assumption (calibration) | Pure equilibrium at 850 degC gives ~0 CH4; real fluidised beds give several vol-% (Puig-Arnavat et al. 2010 discuss such non-equilibrium corrections). Values give 3 / 6 / 8 vol-% dry CH4. |
| A7 | Water-gas shift | at equilibrium at outlet temperature | model choice | Quasi-equilibrium approach, e.g. Zainal et al. (2001). |
| A8 | Higher hydrocarbons, tar, H2S, NH3 | neglected | simplification | Tar/dust removal assumed complete (Woolcock & Brown 2013 review the technologies). |
| A9 | Gasifier heat loss | 2 % of biomass LHV input | assumption | |
| A10 | Steam-to-dry-biomass ratio | oxyfuel 0.5, DFB 0.75 kg/kg; steam at 300 degC | assumption | |
| A11 | O2 purity (oxyfuel) | 99.5 mol-% (balance Ar) | assumption | |
| A12 | DFB combustor heat-transfer efficiency | 0.80 (char + recycled gas LHV -> gasifier) | assumption | |
| A13 | Gas cleaning | cool to 200 degC with 85 % heat recovery, quench to 40 degC, saturation with water, 0.2 bar pressure drop recovered by a blower (eta 0.70) | assumption | |
| A14 | Methanation | equilibrium of CO methanation + WGS at 300 degC, 10 bar (last stage); 2-stage intercooled compression eta 0.75; no carbon deposition; no CO2 removal | assumption / model choice | Fixed-bed methanation conditions reviewed by Rönsch et al. (2016), Kopyscinski et al. (2010). |
| A15 | H2 make-up | partial = 50 %, full = 100 % of the stoichiometric make-up to module M = (H2-CO2)/(CO+CO2) = 3 | definition | |
| A16 | Electrolyser specific energy | 55 kWh/kg H2 (system) | assumption | Inside the ranges reviewed by Buttler & Spliethoff (2018) and Schmidt et al. (2017). |
| A17 | Electrolyser O2 | stoichiometric co-product, 7.94 kg O2/kg H2 | stoichiometry | |
| A18 | ASU specific energy | 0.25 kWh/kg O2 | assumption (order of magnitude) | To be replaced by a referenced value in the scenarios repo. |
| A19 | Burner | complete combustion, lambda = 1.10, combustion air and fuel at 25 / 40 degC, no dissociation | assumption | |
| A20 | Useful heat | flue gas leaves the process at the process temperature (no recuperation) | assumption | Gives the "available heat" of the fuel at that temperature. |
| A21 | Natural-gas reference | 93 % CH4, 4 % C2H6, 1 % C3H8, 1 % CO2, 1 % N2 | assumption | Illustrative H-gas composition, not a measured analysis. |
| A22 | Thermodynamic data | NASA 7-coefficient polynomials | data | GRI-Mech 3.0 thermo data (as shipped with Cantera `gri30.yaml`). |
| A23 | Water saturation pressure | Buck equation (1996 revision) | correlation | |

## Reference design point

10 MW of process heat at a flue-gas exit (process) temperature of 1000 degC, steady state.
