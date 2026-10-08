import cantera as ct
import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar
from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt
import math

# def quadratic_solver(x0,y0,m0,x1,y1,m1):

#     A = np.array([[x0**2,x0,1],
#                  [x1**2,x1,1],
#                  [2*x0,1,0],
#                  [2*x1,1,0]])

#     B = np.array([y0,y1,m0,m1])

#     a,b,c = np.linalg.lstsq(A,B, rcond=None)[0]

#     return a, b, c

def mass_flux(gas: ct.Solution, pressure_ratio: float):

    gas1 = ct.Solution("Combustion.yaml")

    gas1.TPX = gas.TPX

    entropy = gas1.entropy_mass
    ho = gas1.enthalpy_mass
    pressure = (gas1.P)*pressure_ratio

    gas1.SP = entropy, pressure

    gas1.equilibrate('SP')

    h = gas1.enthalpy_mass

    v = (2*(ho - h))**0.5

    density = gas1.density

    mass_flux = density*v

    return mass_flux

# def Nozzle_Inlet_Attributes(gas: ct.Solution, pressure_ratio: float):

#     gas1 = ct.Solution("Combustion.yaml")

#     gas1.TPX = gas.TPX

#     ho = gas1.enthalpy_mass

#     pressure = gas1.P*pressure_ratio

#     gas1.P = pressure

#     gas1.equilibrate('SP')

#     h = gas1.enthalpy_mass
#     cp_mass = gas1.cp_mass
#     cv_mass = gas1.cv_mass
#     gamma = cp_mass/cv_mass
#     T = gas1.T
#     density = gas1.density
#     velocity = (2*(ho - h))**0.5
#     thermal_conductivity = gas1.thermal_conductivity

#     return gamma, T, density, velocity, thermal_conductivity

 

def mach_number(gas: ct.Solution, pressure_ratio: float):

    gas1 = ct.Solution("Combustion.yaml")
    
    gas1.TPX = gas.TPX

    entropy = gas1.entropy_mass
    pressure = gas1.P*pressure_ratio
    ho = gas1.enthalpy_mass
    gas1.SP = entropy , pressure

    gas1.equilibrate('SP')

    h = gas1.enthalpy_mass

    T = gas1.T

    gamma = gas1.cp_mass / gas1.cv_mass

    R = ct.gas_constant/gas1.mean_molecular_weight

    v = (2*(abs(ho-h)))**0.5

    if pressure_ratio ==1:

        v = 0

    mach_number = v/((gamma*R*T)**0.5)

    return mach_number


# def SMD_Plain_Jet (surface_tension, gas_density, lox_tube_diameter, liquid_viscosity, liquid_density, 
#          mass_flow_gas, mass_flow_liquid, gas_velocity, liquid_velocity):

#     GLR = mass_flow_gas/mass_flow_liquid

#     mu_l = liquid_viscosity

#     d = lox_tube_diameter

#     rho_g = gas_density

#     rho_l = liquid_density

#     u_g = gas_velocity

#     u_l = liquid_velocity

#     sigma = surface_tension

#     u_R_squared = (u_g - u_l)**2

#     part_1_A = (sigma/(rho_g*u_R_squared*d))**0.4

#     part_1_B = (1 + ((GLR)**-1))**0.4

#     part_1 = 0.48*d*part_1_A*part_1_B

#     part_2_A = ((mu_l**2)/(sigma*rho_l*d))**0.5

#     part_2_B = 1 + ((GLR)**-1)

#     part_2 = 0.15*d*part_2_A*part_2_B

#     SMD = part_1+part_2

#     return SMD

def Heat_of_Vaporization_LOX_Qdot(Pressure: float, mass_flow_oxidizer: float):

    P = Pressure*ct.one_atm

    h_liquid = CP.PropsSI('H','P',P,'Q',0,'O2')

    h_gas = CP.PropsSI('H','P',P,'Q',1,'O2')

    h_vap = h_gas-h_liquid

    h_vap_dot = h_vap*mass_flow_oxidizer

    return h_vap_dot

def Caloric_Heating_Value_Hdot(H_O2:float, H_Ethanol:float,mass_flow_oxidizer:float,mass_flow_fuel:float,    ##J/s
                          total_mass_flow:float, MM_fuel:float, MM_oxidizer:float, mass_fracs:float):   

    gas1 = ct.Solution("Combustion.yaml")

    
    Y = mass_fracs

    Tref = 298.15

    Pref = ct.one_atm

    gas1.TPY = Tref, Pref, Y

    enthalpy_mass_ref = gas1.enthalpy_mass

    CHVR =(total_mass_flow*enthalpy_mass_ref) - ((mass_flow_fuel*H_Ethanol*(1/MM_fuel))+(mass_flow_oxidizer*H_O2*(1/MM_oxidizer)))

    return abs(CHVR)



def Gas_Film(gas:ct.Solution):

    gas_film = ct.Solution("Combustion.yaml")

    pressure = gas.P

    Chamber_Temperature = gas.T

    Ts = CP.PropsSI('T', 'P', pressure, 'Q', 0, 'O2')
    
    Y_prod = gas.Y
    
    Y_full_O2 = np.zeros(gas.n_species)
    
    Y_full_O2[gas.species_index('O2')] = 1
    
    Y_ref = Y_full_O2 + ((Y_prod - Y_full_O2)/3)

    T = Ts+((Chamber_Temperature-Ts)/3)

    gas_film.TPY = T, pressure, Y_ref

    return gas_film


def Spalding_Number(H:float, Q:float, Mo:float, Chamber_Temperature:float, film_gas: ct.Solution):

    gas = ct.Solution("Combustion.yaml")

    gas_1 = ct.Solution("Combustion.yaml")

    Pressure = film_gas.P

    Y_ref = film_gas.Y

    gas_1.TPY = Chamber_Temperature, Pressure, Y_ref

    gas_ref_chamber_enthalpy = gas_1.enthalpy_mass

    O2_Saturation_Temperature = CP.PropsSI('T', 'P', Pressure, 'Q', 0, 'O2')
           
    gas.TPY = O2_Saturation_Temperature, Pressure*ct.one_atm, Y_ref

    enthalpy_ref_sat=gas.enthalpy_mass

    deltaH = gas_ref_chamber_enthalpy - enthalpy_ref_sat


    # r = (coef_fuel*MM_fuel)/(coef_oxidizer*MM_oxidizer)

    B = ((H/Q) * Mo) + (deltaH/Q)

    return B


def Drag_Actuation_S(film_gas: ct.Solution, Spalding_Number:float):

    mu = film_gas.viscosity
    cp = film_gas.cp_mass
    k = film_gas.thermal_conductivity
    B = Spalding_Number

    Pr = cp*mu/k

    S = (9*Pr)/(2*(np.log(1+B)))

    return S


##User Input##
Thrust_kN = float(input("Enter Target Thrust (kN): "))
pressure = float(input("Enter Chamber Pressure (atm): "))
nozzle_exit_pressure = float(input("Enter Nozzle Exit Pressure (atm): "))
ambient_pressure = float(input("Enter Ambient Pressure (atm): "))
Mc = float(input("Enter Target Mach Number at Nozzle Entrance: "))
phi = float(input("Enter fuel-oxidizer equivalence ratio: "))
# Lstar = float(input("Enter Lstar (m): "))
# Lc_Lt = float(input("Enter length of chamber cylinder to length of chamber cone ratio: "))
length_percent = float(input("Enter bell-nozzle length (percentage of standard 15 degree cone length) (%): "))
theta_i = float(input("Enter initial nozzle angle (deg) post inflection point: "))
theta_e = float(input("Enter nozzle exit angle (deg): "))


pi = math.pi
nozzle_refrence_angle = 15
throat_arc_divergence_coef = 0.4
# throat_arc_convergence_coef = 1.5
print("\n")

Thrust=Thrust_kN*1000

chamber_pressure = pressure*ct.one_atm

exit_pressure = nozzle_exit_pressure*ct.one_atm

ambient_pressure = ambient_pressure*ct.one_atm

##Enthalpy of Reaction


gas = ct.Solution("Combustion.yaml")

gas_burned_copy  = ct.Solution("Combustion.yaml")

gas.TP=298.15, ct.one_atm

all_enthalpies = gas.partial_molar_enthalpies

###O2###
Hf_O2_LOX=CP.PropsSI('Hmolar','T', 90, 'P', chamber_pressure, 'O2')

Hf_O2_Gas=CP.PropsSI('Hmolar_idealgas','T', 298.15, 'P', ct.one_atm, 'O2')

deltaH_O2=Hf_O2_LOX-Hf_O2_Gas ##!!!##

O2_index=gas.species_index('O2')

Hf_O2=all_enthalpies[O2_index]/1000 ##!!!##

###Ethanol###

Hf_ethanol_gas= CP.PropsSI('Hmolar_idealgas','T', 298.15, 'P', chamber_pressure, 'Ethanol')

Hf_ethanol_liquid= CP.PropsSI('Hmolar','T', 298.15, 'P', chamber_pressure, 'Ethanol') 

deltaHf_ethanol=Hf_ethanol_liquid-Hf_ethanol_gas ##!!!##

index_ethanol = gas.species_index('C2H5OH')

Hf_ethanol= all_enthalpies[index_ethanol]/1000 ##!!!##


##H2O##

# index_water = gas.species_index('H2O')

# Hf_water= all_enthalpies[index_water]/1000 ##!!!##



# ##CO2##

# index_CO2=gas.species_index('CO2')

# Hf_CO2= all_enthalpies[index_CO2]/1000 ##!!!##

##deltaHrxn##

gas.set_equivalence_ratio(phi, 'C2H5OH', 'O2')

molefrac_C2H5OH = gas.X[gas.species_index('C2H5OH')]
molefrac_O2 = gas.X[gas.species_index('O2')]

MM_ETHANOL = 0.046069 #kg/mol

MM_O2 = 0.031998 #kg/mol

mass_ethanol = molefrac_C2H5OH*MM_ETHANOL   
mass_O2 = molefrac_O2*MM_O2

O2_ethanol_mass_ratio = mass_O2/mass_ethanol

#deltH_rxn=(((2*Hf_CO2)+(3*Hf_water))-((Hf_ethanol)+(3*Hf_O2))) - (deltaHf_ethanol+(3*deltaH_O2))

#print(f"Enthalpy of Reaction at {chamber_pressure/101325} atm: {(deltH_rxn/1000):.2f} KJ/mol for Ethanol\n") #print(deltH_rxn)


##Equilibrate to find Chamber Conditions##
# db = {s.name: s for s in ct.Species.list_from_file('nasa_gas.yaml')}
# gas = ct.Solution(thermo='ideal-gas',
#                 species=[db[n] for n in ['C2H5OH', 'O2', 'CO2', 'H2O', 'H2', 'CO', 'H', 'O', 'OH']])

# phi = 1
# fuel = 'C2H5OH'
# oxidizer = 'O2'

# gas.set_equivalence_ratio(phi, fuel, oxidizer)

gas.equilibrate('TP')

H_O2 = (Hf_O2 + deltaH_O2)
H_ethanol = (Hf_ethanol + deltaHf_ethanol)

H_total = (molefrac_C2H5OH*H_ethanol) + (molefrac_O2*H_O2)

M_total = (mass_ethanol) + (mass_O2)

H_in = H_total/M_total

gas.HP = H_in, chamber_pressure

gas.equilibrate('HP')

# specific_gas_constant_chamber = ct.gas_constant/gas.mean_molecular_weight # needed for c* calculation

# K_chamber = gas.cp/gas.cv

Tc=gas.T

mass_fracs_product = gas.Y

# gas_burned_copy.TPX = gas.TPX

film_gas = Gas_Film(gas)

mass_specific_heat_chamber = gas.cp_mass

cv_specific_heat_chamber = gas.cv_mass

R = ct.gas_constant/gas.mean_molecular_weight

gamma = mass_specific_heat_chamber/cv_specific_heat_chamber

specfic_heat_film = film_gas.cp_mass

thermal_conductivity_film = film_gas.thermal_conductivity

#Pressure ratio and mass-flux at Nozzle inlet#
nozzle_inlet_pressure_ratio = brentq(lambda pressure_ratio: ( mach_number(gas, pressure_ratio) - Mc), 0.3, 1)

nozzle_inlet_mass_flux = mass_flux(gas, nozzle_inlet_pressure_ratio)

index_H2=gas.species_index('H2')
index_CO=gas.species_index('CO')
index_ethanol=gas.species_index('C2H5OH')
index_water=gas.species_index('H2O')
index_CO2=gas.species_index('CO2')
index_O2=gas.species_index('O2')

print(f"Chamber Thermo-Chemical Properties at {chamber_pressure/101325:.3f} atm:\n")
print(f"Chamber Temperature: {Tc:.3f} K")
print(f"Gamma Products: {gamma:.4f}")
print(f"R of Chamber: {R:.4f}")
print(f"Specific Heat of Film: {specfic_heat_film:.4f}")
print(f"Thermal Conductivity of Film: {thermal_conductivity_film:.4f}")
print(f"H2 mole fraction: {gas.X[index_H2]:.4f}")
print(f"CO mole fraction: {gas.X[index_CO]:.4f}")
print(f"CO2 mole fraction: {gas.X[index_CO2]:.4f}")
print(f"C2H5OH mole fraction: {gas.X[index_ethanol]:.4f}")
print(f"H2O mole fraction: {gas.X[index_water]:.4f}")
print(f"OH mole fraction: {gas.X[gas.species_index('OH')]:.4f}")
print(f"H mole fraction at: {gas.X[gas.species_index('H')]:.4f}")
print(f"O mole fraction at: {gas.X[gas.species_index('O')]:.4f}")
print(f"O2 mole fraction at: {gas.X[index_O2]:.4f}\n")


#Nozzle Throat Thermal properties##

mass_flux_throat_pressure = minimize_scalar(lambda pressure_ratio: -mass_flux(gas, pressure_ratio), bounds=(0.3, 0.9))

throat_pressure_ratio = mass_flux_throat_pressure.x

mass_flux_throat = -mass_flux_throat_pressure.fun

cstar = chamber_pressure/mass_flux_throat

print(f"Throat Thermo-Chemical Properties:")
print(f"Throat Pressure ratio: {throat_pressure_ratio:.3f}")
print(f"Throat Mass Flux: {mass_flux_throat:.2f} kg/(m^2*s)")
print(f"C*: {cstar:.3f} m/s\n")


#Nozzle Inlet and Exit Properties

entropy = gas.entropy_mass

enthalpy_chamber = gas.enthalpy_mass

gas.SP = entropy, exit_pressure

gas.equilibrate('SP')

Te = gas.T

enthalpy_exit_nozzle = gas.enthalpy_mass

ue = (2*(enthalpy_chamber - enthalpy_exit_nozzle))**0.5

desnsity_exit_nozzle = gas.density

nozzle_exit_mass_flux = desnsity_exit_nozzle*ue

print(f"Nozzle Exit and Inlet Thermal Properties:\n")
print(f"Exit velocity (ue) = {ue:.3f} m/s")
print(f"Nozzle Exit Temperature (Te) = {Te:.3f} K")
print(f"Nozzle Exit Mass Flux = {nozzle_exit_mass_flux:.3f} kg/(m^2*s)")
print(f"Nozzle Inlet Pressure Ratio: {nozzle_inlet_pressure_ratio:.3f}\n")

##Mass Flow and Throat, Nozzle Exit, Chamber Area Calculation##

Pc_cstar = chamber_pressure/cstar

mft_mfn = mass_flux_throat/nozzle_exit_mass_flux

mft_mfc = mass_flux_throat/nozzle_inlet_mass_flux

At = Thrust/((Pc_cstar*ue) + ((exit_pressure-ambient_pressure)*mft_mfn))

mdot = (chamber_pressure*At)/cstar

mdot_ethanol = mdot/(O2_ethanol_mass_ratio+1)

mdot_O2 = mdot - mdot_ethanol

Ac = At*mft_mfc

mft_mfe = mass_flux_throat/nozzle_exit_mass_flux

Ae = At*mft_mfe

nozzle_exit_radius = (Ae/pi)**0.5

chamber_radius = (Ac/pi)**0.5

Throat_radius = (At/pi)**0.5

Ae_At = Ae/At

chamber_mass_flux = mdot/Ac

CHV_Ethanol = (Caloric_Heating_Value_Hdot(H_O2, H_ethanol, mdot_O2, mdot_ethanol, mdot, MM_ETHANOL,
                                         MM_O2, mass_fracs_product))/mdot_ethanol

Q_O2 = (Heat_of_Vaporization_LOX_Qdot(pressure, mdot_O2)/mdot_O2)

SpaldingNumber = Spalding_Number(CHV_Ethanol,Q_O2,1,Tc, film_gas)

S = Drag_Actuation_S(film_gas,SpaldingNumber)

print("Specifications:\n")
print(f"Total Mass Flow Rate (kg/s) = {mdot:.3f}")
print(f"Mass Flow Rate Ethanol (kg/s): {mdot_ethanol:.3f}")
print(f"Mass Flow Rate O2 (kg/s): {mdot_O2:.3f}\n")
print(f"Throat Area (m^2): {At:.7f}")
print(f"Throat Radius (m): {Throat_radius:7f}\n")
print(f"Chamber Cross Section Area (m^2): {Ac:7f}")
print(f"Chamber Radius (m): {chamber_radius:7f}\n")
print(f"Nozzle Exit Area (m^2): {Ae:7f}")
print(f"Nozzle Exit Radius (m): {nozzle_exit_radius:7f}\n")
print(f"Nozzle Exit to Throat Aspect Ratio: {Ae_At:7f}\n")
print("Spalding Number and Drag Acutation Coefficient:\n")
print(f"Spalding Number: {SpaldingNumber:.4f}")
print(f"Drag Actuation Coefficient: {S:.4f}")


#Chamber Volume Calculation

# Vc = At*Lstar

# print(f"Chamber Volume (m^3): {Vc:7f}")

# #Chamber Length Calculation

# At_Ac = At/Ac

# Lt = Vc/(Ac*((Lc_Lt+((1+((At_Ac)**0.5)+At_Ac)/3))))

# Lc = Lc_Lt*Lt

# print(f"Chamber Cone Length (m): {Lt:.6f}")
# print(f"Chamber Cylinder Length (m): {Lc:.6f}\n")

# input("Once you have determined the proper specs of the nozzle based on the thorat to nozzle aspect ratio\nPress Enter to continue...")

#Nozzle Length Caluclation

L_cone = (nozzle_exit_radius-Throat_radius)/math.tan((2*pi/360)*nozzle_refrence_angle)

L_bell_nozzle = L_cone * (length_percent/100)

theta_i = theta_i * (2*pi/360)

theta_e = theta_e * (2*pi/360)

xi = (throat_arc_divergence_coef * Throat_radius)*math.sin(theta_i)

yi = Throat_radius +(throat_arc_divergence_coef*Throat_radius*(1-math.cos(theta_i)))

m_i = math.tan(theta_i)

xe = (L_bell_nozzle)

ye = nozzle_exit_radius 

m_e = math.tan(theta_e)

# a,b,c = quadratic_solver(xi,yi,m_i,xe,ye,m_e)

# x = np.linspace(xi,xe,100)

# y = (a*(x**2)) + (b*x) + c 


Px = (yi-ye+(m_e*xe)-(m_i*xi))/(m_e-m_i)

Py = (m_i*(Px-xi)) + yi

t = np.linspace(0,1,200)

x = (((1-t)**2)*xi) + ((2*(1-t)*t)*Px) + ((t**2)*xe)

y = (((1-t)**2)*yi) + ((2*(1-t)*t)*Py) + ((t**2)*ye)

plt.figure(1)
plt.plot(x,y)
# plt.plot(xe,ye)
ann1 = plt.annotate(f'End point  ({xe*1000:.4f}, {ye*1000:.4f}) mm',
                 (xe, ye), textcoords='offset points', xytext=(6, 6), fontsize=8)
# plt.plot(xi,yi)
ann2 = plt.annotate(f'Start point  ({xi*1000:.4f}, {yi*1000:.4f}) mm',
                 (xi, yi), textcoords='offset points', xytext=(6, 6), fontsize=8)
# plt.plot(Px,Py)
ann3 = plt.annotate(f'Control point  ({Px*1000:.4f}, {Py*1000:.4f}) mm',
                 (Px, Py), textcoords='offset points', xytext=(6, 6), fontsize=8)
ann1.draggable()
ann2.draggable()
ann3.draggable()
plt.title('Nozzle Geometry')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.axis('equal')
plt.grid(True)
plt.show()














