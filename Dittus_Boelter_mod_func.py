# -*- coding: utf-8 -*-
'''
hconv_in = Trnsys534(L,rhof,muf,Prf,kf,betaf,nuf,alphaf)

** Inputs
Di: [m] Diametre interne du diametre
DHX: [m] Diametre de l'echangeur de type coil
m_in: [kg/s] Debit dans l'echangeur
pitch: [m] Distance entre deux "loop" du coil exchanger
Cp : [J/(kg K)] Chaleur massique
kf: [W/(m K)] Thermal conductivity
muf: [Ns/m2] Fluid Dynamic viscosity
Prf: [-] Prandtl number

** Output
hconv : [W/m^2] Convection heat transfer coefficient
'''
import numpy as np
def Dittus_Boelter_mod(Di,DHX,pitch,m_in,Cp,kf,muf,Prf):

    ReHX = 4*m_in/(np.pi*Di*muf);
    PrHX= Cp*muf/kf;
    Recrit = 20000*(Di/DHX)**0.32;

    # Cas laminaire
    if ReHX < Recrit:
        HE = ReHX*(Di/DHX)**0.5/(1+(pitch/(np.pi*Di))**2);
        NuHX = ((48/11+(51/11)/(1+1342/PrHX/HE**2))**3+1.816*(HE/(1+1.15/PrHX))**1.5)**(1/3);
    # Cas turbulent
    else:
        NuHX = 0.023*ReHX**0.85*PrHX**0.4*(Di/DHX)**0.1;

    hconv_in = NuHX*kf/Di;    # [W/m^2] Convection heat transfer coefficient

    return hconv_in


