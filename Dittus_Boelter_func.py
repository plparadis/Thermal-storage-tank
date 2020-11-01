# -*- coding: utf-8 -*-
'''
hconv_in = Dittus_Boelter(Vz,D,P,T,fluid)
Convection heat transfer calculation based on Dittus-Boelter correlation

** Inputs
Vz: [m/s] Fluid velocity
D: [m] Pipe hydraulic diameter
P: [kPa] Pressure of fluid
Tk: [K] Temperature of the fluid
fluid: [-] heat transfert fluid

** Output
hconv_in: [W/m^2] Convection heat transfer coefficient
---Properties--------
rho: [kg/m3] Density
mu: [Ns/m2] Dynamic viscosity
Pr: [-] Prandtl number
k: [W/(m K)] Thermal conductivity

'''
def Dittus_Boelter(Vz,D,rho,mu,Pr,k,HeatTransfer):

    ReD = rho*Vz*D/mu;   # [-] Reynolds number

    if HeatTransfer == 'Heating':
        n = 0.4;
    else: # HeatTransfer == 'Cooling'
        n = 0.3;

    NuD = 0.023*ReD**(4/5)*Pr**n; # [-] Nusselt number

    hconv_in = NuD*k/D;    # [W/m^2] Convection heat transfer coefficient
    return hconv_in








