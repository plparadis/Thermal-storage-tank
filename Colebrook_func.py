# -*- coding: utf-8 -*-
'''
dpdz_friction = Colebrook(Vz,D,rugosity_ratio,rhof,muf)
Pressure drop calculation du to friction based on Colebrook correlation

** Inputs
Vz: [m/s] Fluid velocity
D: [m] Pipe inside diameter
rugosity_ratio: [-] ration of absolute rugosity to inside diameter of pipe
-------Properties----------
rhof: [kg/m3] Fluid Density
muf: [Ns/m2] Fluid Dynamic viscosity

** Output
dpdz_friction: [Pa/m] Frictional Pressure drop

'''
import numpy as np
import sys
from scipy.optimize import bisect

eps = sys.float_info.epsilon


def Colebrook(Vz, D, rugosity_ratio, rhof, muf):
    ReD = np.multiply(rhof, Vz) * D / muf  # [-] Reynolds number
    f = np.zeros([np.size(ReD), 1])  # need to initialize that vector

    f[ReD < 2300] = 64 / ReD[ReD < 2300]  # si l'écoulement est laminaire

    fun = lambda x, ReD: -1 / np.sqrt(x) - 2 * np.log10(
        (2.51 / (ReD * np.sqrt(x))) + (rugosity_ratio / 3.7))  # on utilise la corrélation de Colebrook
    f[ReD >= 2300] = np.array(list(map(lambda a: bisect(lambda x: fun(x, a), eps, 0.1), ReD[ReD >= 2300])))
    dpdz_friction = np.multiply(np.multiply(f, rhof), Vz) ** 2 / (2 * D)  # [Pa/m] liquid only pressure drop
    return dpdz_friction
