#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Atmosphere test - tests atmospheres.py
Created on Wed Nov 30 17:54:20 2016

@author: tr1010
"""
import numpy as np
import atmospheres as atmo
import matplotlib.pyplot as plt

rho = np.zeros((8,))
P = np.zeros((8,))
T = np.zeros((8,))
mfp = np.zeros((8,))
eta = np.zeros((8,))
MolW = np.zeros((8,))

heights = np.array([100, 110, 120, 150, 160, 170, 190, 230])*1e3 + 6378.137e3
                   
for i in range(0,8):
    rho[i], P[i], T[i], mfp[i], eta[i], MolW[i] = atmo.US62_76(heights[i])


plt.figure()
plt.semilogx(rho,heights)
plt.figure()
plt.plot(T,heights)