#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Atmospheres: functions which calcualte atmospheric quantities using different
models

Created on Tue Nov 29 11:45:15 2016

@author: tr1010
"""
import numpy as np

# converts geometric altitude to geopotential altitude
def me2po(E,Z):
    H = E*Z/(E + Z)
    return H
    
# US standard Atmosphere
def US62_76(r):
    #Some constants:
    E = 6378.137e3 #6356.0
    Na = np.float64(6.0220978e23)
    sig = np.float64(3.65e-10)
    
    # Sea level standard values:
    P0 = 101325.0 #Pa
    T0 = 288.15 #K
    M = 28.9644 # g/mol
    R = 8.31432  # J/mol-K
    g0 = 9.806658 # m/s2
    GM_R = g0*M/R # GM/R K/km
    Z = (r - E)*1e-3 # convert radius in m to altitude in km
    H = me2po(E,Z) # geopotential altitude
    
    BLdT = np.array([[0., 11., 20., 32., 47., 51., 71., me2po(E,86.), 
                      me2po(E,100.), me2po(E,110.), me2po(E,120.), me2po(E,150.), 
                      me2po(E,160.), me2po(E,170.), me2po(E,190.), me2po(E,230.), 
                      me2po(E,300.), me2po(E,400.), me2po(E,500.), me2po(E,600.), 
                      me2po(E,700.)],
                     [0., -6.5, 0., 1., 2.8, 0., -2.8, -2., 1.693, 5., 10., 20.,
                      15., 10., 7., 5., 4., 3.3, 2.6, 1.7, 1.1]])
    BLT = np.zeros((21,))
    BLP = np.zeros((21,))
    BLT[0] = T0
    BLP[0] = P0 

    for i in range(0, 20):
        # Calculate base temperatures
        BLT[i+1] = BLT[i] + BLdT[1,i+1]*(BLdT[0,i+1] - BLdT[0,i])
        
        # Calculate base pressures
        if (i+1 == 0) or (i+1 == 2) or (i+1 == 5):
            BLP[i+1] = BLP[i]*np.exp(-GM_R*(BLdT[0,i+1] - BLdT[0,i])/BLT[i])
        else:
            BLP[i+1] = BLP[i]*((BLT[i] + BLdT[1,i+1]*(BLdT[0,i+1] - BLdT[0,i]))/BLT[i])**(-GM_R/BLdT[1,i+1])
        
        # Calculate values at requested altitude
        if H > BLdT[0,i] and H <= BLdT[0,i+1]:
            # Temp
            T = np.float64(BLT[i] + BLdT[1,i+1]*(H - BLdT[0,i]))
            
            # Press
            if i+1 == 0 or i+1 == 2 or i+1 == 5:
                P = np.float64(BLP[i]*np.exp(-GM_R*(H - BLdT[0,i])/BLT[i]))
            else:
                P = np.float64(BLP[i]*((BLT[i] + BLdT[1,i+1]*(H - BLdT[0,i]))/BLT[i])**(-GM_R/BLdT[1,i+1])) 
            
            # Density
            rho = np.float64(M*1e-3*P/(R*T))
            mfp = np.float64(M/(2**0.5*np.pi*sig**2*rho*Na)) # mean free path
            eta = np.float64(1.458e-6*T**1.5/(T + 110.4)) # dynamic viscosity via sutherland law

    return rho, P, T, mfp, eta