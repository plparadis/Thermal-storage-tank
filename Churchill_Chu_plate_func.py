# -*- coding: utf-8 -*-
"""
hconv_in = Churchill_Chu2(L,rhof,muf,Prf,kf,betaf,nuf,alphaf)
Natural Convection heat transfer calculation for a Vertical Plate based on Churchill and Chu correlation

** Inputs
L: [m] Height of the Vertical Plate
rhof: [kg/m3] Fluid Density
muf: [Ns/m2] Fluid Dynamic viscosity
Prf: [-] Prandtl number
kf: [W/(m K)] Thermal conductivity
betaf: [1/K]  Volumetric expansivity (beta)
alphaf: [m^2/s] Thermal diffusivity

**Mode :
  1 = Free convection for Vertical plate (9.26)
  2 = Free convection for an upper Horizontal plate
  3 = Free convection for a lower Horizontal plate

**Output
hconv : [W/m^2] Convection heat transfer coefficient
%-------------------------------------------------------------------------%
g = 9.81; % [m/s^2] gravitational acceleration
%--------------------------------------------------------------------------
"""


def Churchill_Chu_plate(L, rhof, Prf, kf, betaf, alphaf, muf, Ts, Tinf, mode):
    g = 9.81  # [m/s^2] gravitational acceleration
    RaL = max(g * betaf * rhof * abs(Ts - Tinf) * L ** 3 / (muf * alphaf), 1e3)  # [-] Rayleigh number (9.25)

    if mode == 1:
        NuL = (0.825 + 0.387 * RaL ** (1 / 6) / (
                    (1 + (0.492 / Prf) ** (9 / 16)) ** (8 / 27))) ** 2  # [-] Nusselt number
    elif mode == 2:
        if Ts > Tinf:
            if RaL > 10 ** 7:
                NuL = 0.15 * RaL ** (1 / 3)
            else:
                NuL = 0.54 * RaL ** (1 / 4)
        else:
            NuL = 0.52 * RaL ** (1 / 5)
    else:  # mode == 3
        if Ts < Tinf:
            if RaL > 10 ** 7:
                NuL = 0.15 * RaL ** (1 / 3)
            else:
                NuL = 0.54 * RaL ** (1 / 4)
        else:
            NuL = 0.52 * RaL ** (1 / 5)

    hconv = NuL * kf / L  # [W/(m^2 K)] Convection heat transfer coefficient

    return hconv
