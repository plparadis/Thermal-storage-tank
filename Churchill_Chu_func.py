# -*- coding: utf-8 -*-
"""
Natural Convection heat transfer calculation based on Churchill and Chu correlation
"""

def Churchill_Chu(D, rhof, Prf, kf, betaf, alphaf, muf, Ts, Tinf):
    """
    Natural Convection heat transfer calculation based on Churchill and Chu correlation
    :param D: [m] Pipe inside diameter
    :param rhof: [kg/m3] Fluid Density
    :param Prf: [-] Prandtl number
    :param kf: [W/(m K)] Thermal conductivity
    :param betaf: [1/K]  Volumetric expansivity (beta)
    :param alphaf: [m^2/s] Thermal diffusivity
    :param muf: [Ns/m2] Fluid Dynamic viscosity
    :param Ts: [°C] Surface temperature
    :param Tinf: [°C] Fluid temperature
    :return hconv_out: [W/m^2] Convection heat transfer coefficient
    """
    g = 9.81  # [m/s^2] gravitational acceleration

    RaD = max(g * betaf * rhof * abs(Ts - Tinf) * D ** 3 / (muf * alphaf), 1000)  # [-] Rayleigh number
    NuD = (0.60 + 0.387 * RaD ** (1 / 6) / ((1 + (0.559 / Prf) ** (9 / 16)) ** (8 / 27))) ** 2  # [-] Nusselt number
    hconv_out = NuD * kf / D  # [W/m^2] Convection heat transfer coefficient
    return hconv_out
