#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Atmospheres: functions which calcualte atmospheric quantities using different
models

Created on Tue Nov 29 11:45:15 2016

@author: tr1010
"""
import sys  
sys.path.append('NRLMSISE-00/Python-NRLMSISE-00-master')
from nrlmsise_00_header import *
from nrlmsise_00 import *
import numpy as np
#import 
#import atmospy76 as atmospy
# converts geometric altitude to geopotential altitude
def me2po(E,Z):
    H = E*Z/(E + Z)
    return H
    
def nrlmsise00(doy,year,sec,alt,g_lat,g_long,lst,f107A,f107,ap):
    
    output = [nrlmsise_output() for _ in range(17)]
    Input = [nrlmsise_input() for _ in range(17)]
    flags = nrlmsise_flags()
    aph = ap_array() # For more detailed ap data (i.e more than daily)
    
    flags.switches[0] = 1
    for i in range(1,24):
        flags.switches[i]=1

    Input.doy=doy
    Input.year=year
    Input.sec=sec
    Input.alt=alt
    Input.g_lat=g_lat
    Input.g_long=g_long
    Input.lst=lst
    Input.f107A=f107A
    Input.f107=f107
    Input.ap=ap
    
    if alt > 500e3:
        gtd7d(Input, flags, output)
    else:
        gtd7(Input, flags, output)
    
    """
    OUTPUT VARIABLES:
    d[0] - HE NUMBER DENSITY(CM-3)
    d[1] - O NUMBER DENSITY(CM-3)
    d[2] - N2 NUMBER DENSITY(CM-3)
    d[3] - O2 NUMBER DENSITY(CM-3)
    d[4] - AR NUMBER DENSITY(CM-3)                       
    d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
    d[6] - H NUMBER DENSITY(CM-3)
    d[7] - N NUMBER DENSITY(CM-3)
    d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
    t[0] - EXOSPHERIC TEMPERATURE
    t[1] - TEMPERATURE AT ALT
    """
    #Now process output to get required values
    kb = 1.38064852e-23 # Boltzmann constant
    Na = 6.022140857e23 # avogadro number
    R0 = kb * Na # universal gas constant
    
    #Molecular weights of different components (kg/kmole)
    molecular_weights = np.zeros(8)
    molecular_weights[0] = 4.002602 #He
    molecular_weights[1] = 15.9994 #O
    molecular_weights[2] = 28.0134 #N2
    molecular_weights[3] = 31.9988 #O2
    molecular_weights[4] = 39.948 #AR
    molecular_weights[5] = 1.00794 #H
    molecular_weights[6] = 14.0067 #N
    molecular_weights[7] = 15.9994 #anomalous O
    
    # Calculate partial pressures
    partial_p = np.zeros(8)
    partial_p[0] = d[0]*k*t[1] #He
    partial_p[1] = d[1]*k*t[1] #O
    partial_p[2] = d[2]*k*t[1] #N2
    partial_p[3] = d[3]*k*t[1] #O2
    partial_p[4] = d[4]*k*t[1] #AR
    partial_p[5] = d[6]*k*t[1] #H
    partial_p[6] = d[7]*k*t[1] #N
    partial_p[7] = d[8]*k*t[1] #anomalous O

    #Assuming perfect gas, calculate atmospheric pressure
    pressure_mixture = np.sum(partial_p)
    
    temperature = t[1]    

    mole_fraction = np.divide(partial_p,P)
    
    molecular_weight_mixture = np.sum(np.multiply(mole_fraction,molecular_weights)) #g/mol
    
    mass_fractions = np.multiply(mole_fraction,
                                 np.divide(molecular_weights,molecular_weight_mixture))
    
    specific_gas_constants = R0/molecular_weights
    
    R_mixture = np.sum(np.multiply(mass_fractions,specific_gas_constants))
    
    number_density_mixture = np.sum(d) - d[5] 
    
    
    
#    mfp = np.float64(MolW*1e-3/(2**0.5*np.pi*sig**2*rho*Na)) # mean free path
#    eta = np.float64(1.458e-6*T**1.5/(T + 110.4)) # dynamic viscosity via sutherland law
#    SoS = np.float64(np.sqrt(1.4*287.085*T))
    
    return rho, P, T, R, mfp, eta, MolW, SoS
    
# US mutant Atmosphere
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