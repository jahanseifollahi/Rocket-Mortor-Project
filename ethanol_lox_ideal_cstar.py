import cantera as ct
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import minimize_scalar 



def ideal_gas():
    
    gas = ct.Solution("Combustion.yaml")
    return gas

def equilibrate(phi, Pressure):

    gas = ideal_gas()

    gas.TP = 298.15, Pressure #gas.TP = 298.15, ct.one_atm

    gas.set_equivalence_ratio(phi, 'C2H5OH', 'O2')

    index_O2 = gas.species_index('O2')
    index_ethanol = gas.species_index('C2H5OH')
    
    O2_mols_frac = gas.X[index_O2]
    ethanol_mols_frac = gas.X[index_ethanol]

    Hin = enthalpy_of_mixture(ethanol_mols_frac, O2_mols_frac, Pressure)

    gas.equilibrate('TP')

    gas.HP = (Hin, Pressure)

    gas.equilibrate('HP')


    # index_O2 = gas.species_index('O2')
    # index_ethanol = gas.species_index('C2H5OH')


    # O2_mols_frac = gas.X[index_O2]
    # ethanol_mols_frac = gas.X[index_ethanol]

    # Hin = enthalpy_of_mixture(ethanol_mols_frac, O2_mols_frac, Pressure)

    # molsC= 2*ethanol_mols_frac
    # molsH=6*ethanol_mols_frac
    # molsO=ethanol_mols_frac+(2*O2_mols_frac)
    # total_mols = molsC + molsH + molsO

    # mol_frac_C = molsC/total_mols
    # mol_frac_H = molsH/total_mols
    # mol_frac_O = molsO/total_mols

    # mols_product_C_lean = 2
    # mols_product_H_lean = 6

    # mols_O_per_O2 = 2
    # mols_O_per_water = 1
    # mols_O_per_CO2 = 2
    # mols_C_per_CO2 = 1
    # mols_H_per_H2O = 2

    # mols_C_per_ethanol = 2
    # mols_H_per_ethanol = 6
    # mols_O_per_ethanol = 1

    # mols_product_C_noEthanol = 2*(mols_C_per_CO2)
    # mols_product_H_noEthanol = 3*(mols_H_per_H2O)
    # mols_product_O_noEthanol = 2*(mols_O_per_CO2)+3*(mols_O_per_water)
    # mols_total_noEthanol = mols_product_C_noEthanol+mols_product_H_noEthanol+mols_product_O_noEthanol
    # mols_total_elements_in_Ethanol = mols_C_per_ethanol+mols_H_per_ethanol+mols_O_per_ethanol




    # if phi<1.0:
    #     mols_O2= (((mols_product_H_lean+mols_product_C_lean+(2*mols_O_per_CO2)+(3*mols_O_per_water))*mol_frac_O)-((2*mols_O_per_CO2)+(3*mols_O_per_water)))/(mols_O_per_O2*(1-mol_frac_O))
    #     gas = ideal_gas()
    #     gas.X = f'CO2:2, H2O:3, O2: {mols_O2}'
    #     gas.HP = Hin, Pressure
    #     gas.equilibrate('HP')
       
    # elif phi>1.0:
    #     mols_ethanol = ((mols_total_noEthanol*mol_frac_O)-(mols_product_O_noEthanol))/(1-(mol_frac_O*mols_total_elements_in_Ethanol))
    #     gas = ideal_gas()
    #     gas.X = f'CO2:2, H2O:3, C2H5OH: {mols_ethanol}'
    #     gas.HP = Hin, Pressure
    #     gas.equilibrate('HP')
       
    # else:
    #     gas=ideal_gas()
    #     gas.X = f'CO2:2, H2O:3'
    #     gas.HP = Hin, Pressure
    #     gas.equilibrate('HP')
        
    
    maximum_mass_flux = minimize_scalar(lambda ratios: -mass_flux(gas,ratios), bounds=(0.1,0.8))

    cstar = Pressure/-(maximum_mass_flux.fun)

    pressure_ratio = maximum_mass_flux.x

    return cstar, pressure_ratio

def mass_flux (gas:ct.Solution, pressure_ratio):

    gas1 = ideal_gas()

    gas1.TPX = gas.TPX

    entropy = gas1.entropy_mass

    pressure = (gas1.P)*pressure_ratio

    h0= gas1.enthalpy_mass

    gas1.SP = entropy, pressure

    gas1.equilibrate('SP')

    h = gas1.enthalpy_mass

    v = (2*(h0 - h))**0.5

    density = gas1.density

    mass_flux = density*v

    return mass_flux



def enthalpy_of_mixture(moles_of_fuel,moles_of_oxidizer,Pressure):

    gas = ideal_gas()

    gas.TP = 298.15, ct.one_atm

    enthalpies = gas.partial_molar_enthalpies

    O2_index = gas.species_index('O2')
    ethanol_index = gas.species_index('C2H5OH')

    H_O2 = (enthalpies[O2_index]/1000) + (CP.PropsSI('Hmolar','T',90,'P',Pressure,'O2') - CP.PropsSI('Hmolar_idealgas','T',298.15,'P',101325,'O2'))

    H_Ethanol = (enthalpies[ethanol_index]/1000) + (CP.PropsSI('Hmolar','T',298.15,'P',Pressure,'Ethanol') - CP.PropsSI('Hmolar_idealgas','T',298.15,'P',101325,'Ethanol'))

    H_total = (moles_of_fuel*H_Ethanol) + (moles_of_oxidizer*H_O2)

    M_total = (moles_of_fuel*0.046069) + (moles_of_oxidizer*(0.031998))

    Hin_per_mass = H_total/M_total

    return Hin_per_mass

phis = np.linspace(0.1,3,200)

pressure = float(input("Enter the pressure in atm: "))

pressure = pressure*ct.one_atm

values =[]
ratios = []

for phi in phis:
    value, ratio = equilibrate(phi,pressure)
    values.append(value)
    ratios.append(ratio)

maximum_value = max(values)
idx_max = values.index(maximum_value)

print(f"The maximum value is {maximum_value:.3f} at phi = {phis[idx_max]:.3f} for pressure ratio {ratios[idx_max]:.3f}")

plt.plot(phis,values)
plt.title('C* vs phi')
plt.xlabel('phi')
plt.ylabel('C*')
plt.grid(True)
plt.show()


    
    





