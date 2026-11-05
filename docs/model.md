# Model description

## Agents and scheduler

`biomass_heat.plant.Plant` is the scheduler / message bus. It owns one instance of each
agent and, for every time step, passes messages in this order:

```
ProcessHeatDemand -> (size chain) -> BiomassFeed -> Gasifier -> GasCleaning -> Methanation
      -> Burner                     Strategy -> Electrolyser -> GasStorage(H2), GasStorage(O2)
                                    AirSeparationUnit (O2 top-up) -> GridMarket (cost)
```

Each agent keeps its own parameters, state (storages: fill level) and a ledger of
cumulative flows. Operating strategies are separate objects implementing
`electrolyser_power(t, plant, h2_need, o2_need, price, h2_storage, o2_storage)`;
`FollowDemand` and a simple `PriceThreshold` are included, optimised strategies follow in
`biomass-heat-electrolysis-dispatch`.

## Gasifier (quasi-equilibrium)

Unknowns CO, CO2, H2, H2O for a given oxidant flow, with CH4 fixed (A6):

* C: CO + CO2 = X_C * C_bio - CH4
* H: 2 H2 + 2 H2O = H_bio + 2 (H2O_moisture + H2O_steam) - 4 CH4
* O: CO + 2 CO2 + H2O = O_bio + H2O_moisture + H2O_steam + 2 O2
* WGS: x_CO2 x_H2 / (x_CO x_H2O) = K_WGS(T)  -> quadratic in CO2

Autothermal modes: O2 (with air N2/Ar for air gasification) is found by bisection such that

H_biomass + H_moisture(l) + H_steam + H_oxidant = H_gas(T_out) + Q_loss

with absolute enthalpies (formation + sensible). Biomass formation enthalpy follows from its
HHV. Char leaves at 25 degC with h = 0 (its chemical energy is lost).

Allothermal (DFB) mode: no oxidant; the deficit is covered by burning the char in the
combustor and, if needed, a fraction of the product gas (efficiency A12).

## Methanation

Two independent reactions, CO + 3 H2 = CH4 + H2O and CO + H2O = CO2 + H2, solved for
equilibrium at (T, p) by nested bisection on the reaction extents. Results were cross-checked
against Cantera's Gibbs minimisation with the same species set during development
(identical mole fractions to 4 decimals).

## Burner and available heat

Adiabatic flame temperature from the enthalpy balance of complete combustion. Available heat
= H_reactants - H_fluegas(T_exit). Fuel demand = Q_process / available heat fraction.

## Key indicators

* Cold gas efficiency (CGE) = LHV of clean syngas / LHV of biomass (as received)
* Fuel conversion efficiency = LHV of fuel gas / (LHV biomass + LHV H2 added)
* Chain heat efficiency = process heat / (LHV biomass + electricity)
* Upper Wobbe index = HHV_vol / sqrt(relative density), dry gas, 0 degC / 1.01325 bar
