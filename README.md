# Ethanol / LOX Rocket Engine Sizing

Python scripts for the preliminary design of a liquid rocket engine burning ethanol with liquid oxygen (LOX). Starting from a target thrust and chamber pressure, they compute the combustion-chamber state, the throat, chamber and nozzle-exit sizes, a bell-nozzle contour, and the injector and characteristic-length (L*) quantities needed to size the chamber.

Combustion is modelled as chemical equilibrium with [Cantera](https://cantera.org/). Real-fluid properties of LOX and ethanol (liquid enthalpy, density, saturation state) come from [CoolProp](http://www.coolprop.org/).

## Files

| File | Type | Purpose |
| --- | --- | --- |
| `ethanol_lox_ideal_cstar.py` | Script | Finds the equivalence ratio that maximises C* at a given chamber pressure |
| `ethanol_lox_engine_spec.py` | Script | Sizes the engine and draws the bell-nozzle contour |
| `Lstar_ChamberVolume.py` | Function library (in progress) | Chamber volume and injector orifice sizing |
| `Combustion.yaml` | Data | Cantera gas model used by all of the scripts |
| `Nozzle_Geometry.png` | Output | Example nozzle contour plot from `ethanol_lox_engine_spec.py` |


## Setup

The scripts need `cantera`, `CoolProp`, `numpy`, `scipy` and `matplotlib`. Run them from the project folder, because `Combustion.yaml` is loaded by relative path:

```
conda activate cantera-env
python ethanol_lox_ideal_cstar.py
python ethanol_lox_engine_spec.py
```

Both scripts ask for their inputs at the prompt.

## Typical workflow

1. Run `ethanol_lox_ideal_cstar.py` to choose an equivalence ratio for your chamber pressure.
2. Run `ethanol_lox_engine_spec.py` with that equivalence ratio to get mass flows, areas, radii and the nozzle contour.
3. Once `Lstar_ChamberVolume.py` is complete, it will use the chamber and film properties from step 2 to calculate the chamber volume and size the injector orifices.

## Script details

### `ethanol_lox_ideal_cstar.py`

Sweeps the fuel-oxidizer equivalence ratio (phi) from 0.1 to 3 in 200 steps at one chamber pressure and reports where the characteristic velocity C* peaks.

**Input:** chamber pressure (atm).

**Output:** the maximum C*, the phi at which it occurs, the throat pressure ratio at that point, and a plot of C* against phi.

For each phi it:

- computes the enthalpy of the incoming propellants, with LOX at 90 K and liquid ethanol at 298.15 K (`enthalpy_of_mixture`);
- equilibrates the mixture at constant enthalpy and pressure to get the adiabatic chamber state (`equilibrate`);
- expands the gas isentropically, in shifting equilibrium, and finds the pressure ratio that maximises mass flux, which is the throat (`mass_flux`);
- sets C* = chamber pressure / throat mass flux.

### `ethanol_lox_engine_spec.py`

The main sizing script. It prints the results in sections and then opens the nozzle contour plot.

**Inputs:**

| Prompt | Units |
| --- | --- |
| Target thrust | kN |
| Chamber pressure | atm |
| Nozzle exit pressure | atm |
| Ambient pressure | atm |
| Target Mach number at the nozzle entrance | – |
| Fuel-oxidizer equivalence ratio | – |
| Bell-nozzle length, as a percentage of a 15° cone | % |
| Initial nozzle angle after the inflection point | deg |
| Nozzle exit angle | deg |

**Outputs:**

- **Chamber properties:** temperature, gamma, gas constant, product mole fractions, and the specific heat and thermal conductivity of the gas film around a LOX droplet.
- **Throat properties:** pressure ratio, mass flux and C*.
- **Nozzle properties:** exit velocity, exit temperature, exit mass flux, and the pressure ratio at the nozzle inlet that gives the target Mach number.
- **Specifications:** total, ethanol and O2 mass flow rates; area and radius of the throat, chamber and nozzle exit; exit-to-throat area ratio.
- **Droplet evaporation terms:** the Spalding number and drag actuation coefficient, used in the L* calculation.
- **Nozzle contour:** a plot of the diverging section, with draggable labels for the start, end and control points in mm.

Functions:

| Function | Returns |
| --- | --- |
| `mass_flux` | Mass flux after isentropic equilibrium expansion to a given pressure ratio |
| `mach_number` | Mach number at a given pressure ratio |
| `Heat_of_Vaporization_LOX_Qdot` | Heat rate needed to vaporise the LOX flow |
| `Caloric_Heating_Value_Hdot` | Heat release rate of the reaction at reference conditions |
| `Gas_Film` | Gas film around a LOX droplet, using the 1/3 rule between the droplet surface and the chamber gas |
| `Spalding_Number` | Spalding transfer number B for the LOX droplet |
| `Drag_Actuation_S` | Drag actuation coefficient S, from the film Prandtl number and B |

The throat area comes from the thrust equation, including the pressure-thrust term. The chamber and exit areas follow from the ratio of throat mass flux to the local mass flux. The contour is a quadratic Bézier curve that leaves the throat arc (radius 0.382 × throat radius) at the initial angle and reaches the exit radius at the exit angle.

### `Lstar_ChamberVolume.py`

This code is still in the process of being completed. It will calculate the volume of the chamber and size the orifices of the injectors. In this model LOX is injected as a liquid jet and ethanol as a saturated vapour.

| Function | Returns |
| --- | --- |
| `pressure_drop` | Injector pressure drop for a chamber pressure |
| `Lox_Injection_Velocity` | LOX injection velocity from that pressure drop |
| `Lox_Injecttion_Properties` | LOX orifice area per element, manifold pressure, LOX density |
| `Ethanol_Injection_Velocity` | Ethanol injection velocity for a momentum-flux ratio J, capped at the speed of sound |
| `Ethanol_Injection_Properties` | Ethanol orifice area per element, manifold pressure, vapour density, injection velocity |
| `SMD_Plain_Jet` | Sauter mean diameter of the LOX spray (plain-jet airblast correlation) |
| `Xi` | Dimensionless droplet evaporation length |
| `Lstar` | Characteristic chamber length L* needed to evaporate the LOX droplets |

Pressures passed to these functions are in atm.

### `Combustion.yaml`

The Cantera input file, converted from the San Diego mechanism. It defines an ideal-gas phase with 12 species: N2, AR, HE, C2H5OH, O2, CO2, H2O, H2, CO, H, O and OH, with mixture-averaged transport. The scripts use it for equilibrium and transport properties only.

## Sources

The L*, Spalding number and drag actuation coefficient calculations (`Lstar`, `Spalding_Number` and `Drag_Actuation_S`) follow:

- [Evaluation of SMD Effects on Characteristic Lengths of Liquid Rocket Engines Using Ethanol/LOx and RP-1/LOx](https://www.researchgate.net/publication/345906474_Evaluation_of_SMD_Effects_on_Characteristic_Lengths_of_Liquid_Rocket_Engines_Using_EthanolLOx_and_RP-1LOx)

The Sauter mean diameter correlation (`SMD_Plain_Jet`) comes from:

- A. H. Lefebvre and D. R. Ballal, *Gas Turbine Combustion: Alternative Fuels and Emissions*, 3rd ed., CRC Press, 2010.
