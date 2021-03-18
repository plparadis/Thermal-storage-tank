# -*- coding: utf-8 -*-
"""
hconv_in = Churchill_Chu(D,rhof,muf,Prf,kf,betaf,nuf,alphaf)
Natural Convection heat transfer calculation based on Churchill and Chu correlation

** Inputs
D: [m] Pipe inside diameter
rhof: [kg/m3] Fluid Density
muf: [Ns/m2] Fluid Dynamic viscosity
Prf: [-] Prandtl number
kf: [W/(m K)] Thermal conductivity
betaf: [1/K]  Volumetric expansivity (beta)
alphaf: [m^2/s] Thermal diffusivity

** Output
hconv_in: [W/m^2] Convection heat transfer coefficient
"""


def Churchill_Chu(D, rhof, Prf, kf, betaf, alphaf, muf, Ts, Tinf):
    g = 9.81;  # [m/s^2] gravitational acceleration

    RaD = max(g * betaf * rhof * abs(Ts - Tinf) * D ** 3 / (muf * alphaf), 1000)  # [-] Rayleigh number
    NuD = (0.60 + 0.387 * RaD ** (1 / 6) / ((1 + (0.559 / Prf) ** (9 / 16)) ** (8 / 27))) ** 2  # [-] Nusselt number
    hconv_out = NuD * kf / D  # [W/m^2] Convection heat transfer coefficient
    return hconv_out
