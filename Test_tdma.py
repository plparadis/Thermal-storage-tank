# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:58:02 2019

@author: plparadis
"""
import numpy as np
from TDMAfunc import TDMA



A = np.array([[1,2,0],[2,2,2],[0,-1,2]])
B = np.array([[2],[12],[5]])


X1 = np.linalg.lstsq(A,B, rcond=None)[0]

a = np.array([[2],[2]])
b = np.array([[2],[-1]])
c = B
d = np.array([[1],[2],[2]])

X2 = TDMA(a,b,c,d)


