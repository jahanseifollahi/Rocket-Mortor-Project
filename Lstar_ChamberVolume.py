import CoolProp.CoolProp as CP
import cantera as ct
import numpy as np


def pressure_drop(Pc):

    delta_P = 80*((10*Pc*ct.one_atm)**0.5)

    return delta_P

def Lox_Injection_Velocity(Pc):

    rho_Lox = CP.PropsSI('D', 'T', 90, 'P', Pc*ct.one_atm, 'O2')

    delta_p = pressure_drop(Pc)

    u_inj = np.sqrt(2*delta_p/rho_Lox)

    return u_inj

def Xi(uinj, ugas, drag_actuation):   ## dimensionless evap length

    Xo = uinj/ugas

    Xi = (Xo + ((3*drag_actuation)/10))/(2+drag_actuation)

    return Xi

def Lox_Injecttion_Properties(discharge_coefficient, uinj, mass_flow, Pc , N):

    rho_Lox = CP.PropsSI('D', 'T', 90, 'P', Pc*ct.one_atm, 'O2')

    A_inj = (mass_flow/N)/(discharge_coefficient*rho_Lox*uinj)

    manifold_pressure = Pc + (pressure_drop(Pc)/ct.one_atm)

    return A_inj , manifold_pressure, rho_Lox

def Ethanol_Injection_Velocity(J, lox_velocity,Pc):

    # gas = ct.Solution("Combustion.yaml")
    
    # gas.X = 'C2H5OH:1'
    
    # sat_T = CP.PropsSI('T', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    # cp = CP.PropsSI('C', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    # cv = CP.PropsSI('O', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    # gamma = cp/cv
    
    # gas.TP = sat_T, Pc*ct.one_atm

    # gamma = gas.cp_mass/gas.cv_mass

    # M = CP.PropsSI('M', 'Ethanol')

    # Ru = CP.PropsSI('gas_constant', 'Ethanol')

    # R = Ru/M

    ethanol_density = CP.PropsSI('D', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    lox_density = CP.PropsSI('D', 'T', 90, 'P', Pc*ct.one_atm, 'O2')

    M_1 = CP.PropsSI('A', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    u_inj = ((J*(lox_density*((lox_velocity)**2)))/(ethanol_density))**0.5

    if u_inj > M_1:

        u_inj = M_1

        print(f"For the given J the injection velocity is grater than Mach 1, the injection velocity has been changed to {u_inj:.3f} m/s")
        
    return u_inj

def Ethanol_Injection_Properties(discharge_coefficient, J, mass_flow, Pc , N):

    # gas = ct.Solution("Combustion.yaml")

    # ethanol_index = gas.species_index('C2H5OH')

    # gas.X = 'C2H5OH:1'

    sat_T = CP.PropsSI('T', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    # gas.TP = sat_T, Pc*ct.one_atm

    ethanol_density = CP.PropsSI('D', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    # rho_Lox = CP.PropsSI('D', 'T', 90, 'P', Pc*ct.one_atm, 'O2')

    Lox_inj_velicity = Lox_Injection_Velocity(Pc)
    
    u_inj = Ethanol_Injection_Velocity(J, Lox_inj_velicity, Pc)

    A_inj = (mass_flow/N)/(discharge_coefficient*ethanol_density*u_inj)

    hf = CP.PropsSI('H', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    s = CP.PropsSI('S', 'P', Pc*ct.one_atm, 'Q', 1, 'Ethanol')

    hin = (0.5*(u_inj**2)) + hf

    P_manifold = CP.PropsSI('P', 'H', hin, 'S', s, 'Ethanol')

    # hf = gas.enthalpy_mass

    # ho = (0.5*(u_inj**2)) - hf

    # gas. = ho

    # manifold_pressure = gas.P

    # delta_p = P_manifold - (Pc*ct.one_atm)

    # delta_p = delta_p/ct.one_atm

    return A_inj, (P_manifold/ct.one_atm), ethanol_density, u_inj

def SMD_Plain_Jet (gas_density, lox_tube_diameter, liquid_viscosity, liquid_density, 
         mass_flow_gas, mass_flow_liquid, gas_velocity, liquid_velocity):

    GLR = mass_flow_gas/mass_flow_liquid

    mu_l = liquid_viscosity

    d = lox_tube_diameter

    rho_g = gas_density

    rho_l = liquid_density

    u_g = gas_velocity

    u_l = liquid_velocity

    sigma = 0.0017783 ## N/m

    u_R_squared = (u_g - u_l)**2

    part_1_A = (sigma/(rho_g*u_R_squared*d))**0.4

    part_1_B = (1 + ((GLR)**-1))**0.4

    part_1 = 0.48*d*part_1_A*part_1_B

    part_2_A = ((mu_l**2)/(sigma*rho_l*d))**0.5

    part_2_B = 1 + ((GLR)**-1)

    part_2 = 0.15*d*part_2_A*part_2_B

    SMD = part_1+part_2

    return SMD

def Lstar(SMD,gamma,product_density,film_specific_heat,film_thermal_conductivity, throat_mass_flux, spalding_number,
          Chamber_Temperature, R, Xi, Pressure): 

    desniy_droplet = CP.PropsSI('D', 'T', 90, 'P', Pressure*ct.one_atm, 'O2')

    ro = SMD/2

    square_root_gammRTc = (gamma*R*Chamber_Temperature)**0.5

    gamma_plus_1 = gamma + 1

    gamma_minus_1 = gamma - 1

    part_1 = Xi*((ro)**2)

    part_2_A = 2/gamma_plus_1

    part_2_B = (throat_mass_flux/(product_density*square_root_gammRTc))**2

    part_2_C = gamma_minus_1/gamma_plus_1

    part_2_D = gamma_plus_1/(2*gamma_minus_1)

    part_2 = (part_2_A+(part_2_B*part_2_C))**part_2_D

    part_3_A = ((film_specific_heat*desniy_droplet)/film_thermal_conductivity)

    part_3_B = ((gamma*R*Chamber_Temperature)**0.5)/np.log(1+spalding_number)

    part_3 = part_3_A*part_3_B

    L_star = part_1*part_2*part_3

    return L_star