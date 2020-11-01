# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:04:05 2019

@author: plparadis
"""

h_convEC[i,0] = Churchill_Chu_plate(paramRes2.Di,propriFluid_EC.rho_inf[i,0],
        propriFluid_EC.Pr_inf[i,0],propriFluid_EC.k_inf[i,0],propriFluid_EC.beta_inf[i,0],
            propriFluid_EC.alpha_inf[i,0],propriFluid_EC.mu_inf[i,0],Tsk_in[i,0],Res2results.Tk[j,1],2); # [W/(m^2) K] Convection heat transfer coefficient in EC
