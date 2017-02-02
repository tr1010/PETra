#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Atmospheres: functions which calcualte atmospheric quantities using different
models

Created on Tue Nov 29 11:45:15 2016

@author: tr1010
"""
import numpy as np
import atmospy76 as atmospy
# converts geometric altitude to geopotential altitude
def me2po(E,Z):
    H = E*Z/(E + Z)
    return H
    
# US standard Atmosphere
def US62_76(r):
    #Some constants:
    E = 6378.137e3 #6356.0e3
    Na = np.float64(6.0220978e23)
    sig = np.float64(3.65e-10)
    
    # Sea level standard values:
    P0 = 101325.0 #Pa
    T0 = 288.15 #K
    M = np.array([28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.9644, 28.962, 28.962,
                    28.88, 28.56, 28.07, 26.92, 26.66, 26.4, 25.85,
                    24.70, 22.66, 19.94, 17.94, 16.84, 16.17]) # Molecular masses with altitude g/mol
    R0 = 8.31432  # J/mol-K
    g0 = 9.806658 # m/s2
    GM_R = g0*M/R0 # GM/R K/km
    Z = (r - E)*1e-3 # convert radius in m to altitude in km
    H = me2po(E,Z) # geopotential altitude

    BLH = np.array([0., 11., 20., 32., 47., 51., 71., me2po(E,86.), 
                    me2po(E,100.), me2po(E,110.), me2po(E,120.), me2po(E,150.), 
                    me2po(E,160.), me2po(E,170.), me2po(E,190.), me2po(E,230.), 
                    me2po(E,300.), me2po(E,400.), me2po(E,500.), me2po(E,600.), 
                    me2po(E,700.)])

    L = np.array([0., -6.5, 0., 1., 2.8, 0., -2.8, -2., 1.693, 5., 10., 20., 15., 
                  10., 7., 5., 4., 3.3, 2.6, 1.7, 1.1])
    BLT = np.zeros((21,))
    BLP = np.zeros((21,))
    BLT[0] = T0
    BLP[0] = P0 

    for i in range(0, 20):
        # Calculate base temperatures
        BLT[i+1] = BLT[i] + L[i+1]*(BLH[i+1] - BLH[i])
        
        # Calculate base pressures
        if (i+1 == 0) or (i+1 == 2) or (i+1 == 5):
            BLP[i+1] = BLP[i]*np.exp(-GM_R[i+1]*(BLH[i+1] - BLH[i])/BLT[i])
        else:
            BLP[i+1] = BLP[i]*((BLT[i] + L[i+1]*(BLH[i+1] - BLH[i]))/BLT[i])**(-GM_R[i+1]/L[i+1])
        
        # Calculate values at requested altitude
        if H > BLH[i] and H <= BLH[i+1]:
            # Molecular weight (interpolate)]
            MolW = M[i] + (M[i+1] - M[i])*(H - BLH[i])/(BLH[i+1] - BLH[i])
            gmrtemp = g0*MolW/R0
            
            # Molecular scale Temperature
            T = np.float64(BLT[i] + L[i+1]*(H - BLH[i]))
            T = MolW*T/M[0] # Convert molecular scale temperature to kinetic temperature

            # Pressure
            if i+1 == 0 or i+1 == 2 or i+1 == 5:
                P = np.float64(BLP[i]*np.exp(-gmrtemp*(H - BLH[i])/BLT[i]))
            else:
                P = np.float64(BLP[i]*((BLT[i] + L[i+1]*(H - BLH[i]))/BLT[i])**(-gmrtemp/L[i+1])) 
            
            # Density
            rho = np.float64(MolW*1e-3*P/(R0*T))
            mfp = np.float64(MolW*1e-3/(2**0.5*np.pi*sig**2*rho*Na)) # mean free path
            eta = np.float64(1.458e-6*T**1.5/(T + 110.4)) # dynamic viscosity via sutherland law
            SoS = np.float64(np.sqrt(1.4*287.085*T))

    return rho, P, T, mfp, eta, MolW, SoS
    
def US76_FORTRAN(alt):
    alt = alt*1e-3
    sig_atmo,del_atmo,the_atmo = atmospy.atmosphere(alt)
    
    #Some constants:
    Na = np.float64(6.0220978e23) # avogadro number
    d = np.float64(3.65e-10)      # molecule diameter
    # Sea level standard values:
    P0 = 101325.0 #Pa
    T0 = 288.15 #K
    R = 287.058
    rho0 = P0/(R*T0)
    asound0 = np.sqrt(1.4*R*T0)
    
    T=T0*the_atmo
    P=P0*del_atmo
    rho=rho0*sig_atmo
    asound=asound0*np.sqrt(the_atmo)
    eta = np.float64(1.458e-6*T**1.5/(T + 110.4)) # dynamic viscosity via sutherland law
    mfp = 8.3145*T /(2**0.5*np.pi*d**2*Na*P)
    
    
    return rho, P, T, mfp, eta, asound