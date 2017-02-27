#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
atmospheres.py includes functions to calculate atmospheric quantities.

Created on Tue Nov 29 11:45:15 2016

@author: tr1010 (Thomas Rees)
"""
import sys  
sys.path.append('atmosphere_models/Python-NRLMSISE-00-master')
from nrlmsise_00_header import *
from nrlmsise_00 import *
import numpy as np
    
def nrlmsise00(doy,year,sec,alt,g_lat,g_long,lst,f107A,f107,ap):
    """
    nrlmsise00 calculates atmospheric quantities using the NRLMSISE-00
    atmosphere published in 2001 by Mike Picone, Alan Hedin, and Doug Drob.
    Originally written in FORTRAN, it was later implemented in C by Dominik 
    Brodowski.
    
    This function calls a Python port of Brodowski's C implementation originally
    written by Joshua Milas in 2013. This software was released under an MIT 
    license (see the license file in the atmosphere_models directory).
    
    The NRLMSISE-00 model uses a number of switches (contained in the flags 
    class) to modify the model output.  At the moment, these defaults are hard-
    wired into PETra. Later revisions will give the user the ability to select
    these switches. For more detailed information about the inputs/outputs/switches
    used in this model, the user is directed to the docstrings of the funcitons
    contained in the model files (norlmsise_00_header.py and nrlmsise_00.py).
    
    Inputs: 
        doy: day of year
        year: year (currently ignored)
        sec: seconds in day
        alt: altitude
        g_lat: geodetic latitude
        g_long: geodetic longitude
        lst: local apparent solar time (hours)
        f107A: 81 day average of F10.7 flux (centred on doy)
        f107: daily f10.7 flux (for previous day)
        ap: magnetic index (daily)
     
    Outputs:
        rho: density at the requested altitude
        pressure_mixture: pressure at the requested altitude
        temperature: temperature at the requested altitude
        R_mixture: the gas constant of the mixture
        mean_free_path: mean free path of the air at the requested altitude. 
                        In contrast to the other outputs of this function, the
                        mean free path calculation assumes a single molecule 
                        gas (assumed to be an 'average' air molecule)
        eta: viscosity (calcualted using Sutherland's law)
        molecular_weight_mixture: the molecular weight of the air at the
                                  requested altitude
        SoS: speed of sound (assume ratio of specific heats is constant 1.4
             everywhere in the atmosphere)
    
    """
    
    output = nrlmsise_output()
    Input = nrlmsise_input()
#    output = [nrlmsise_output() for _ in range(17)]
#    Input = [nrlmsise_input() for _ in range(17)]
    flags = nrlmsise_flags()
    aph = ap_array() # For more detailed ap data (i.e more than daily)
    
    flags.switches[0] = 1 # to have results in m rather than cm
    for i in range(1,24):
        flags.switches[i]=1

    # below 80 km solar & magnetic effects not well established so set to defaults
    if alt < 80e3:
        f107 = 150.
        f107A = 150.
        ap = 4.
    
    # fill out Input class
    Input.year=year
    Input.doy=doy
    Input.sec=sec
    Input.alt=alt*1e-3 #change input to km
    Input.g_lat=g_lat*180/np.pi
    Input.g_long=g_long*180/np.pi
    Input.lst=lst
    Input.f107A=f107A
    Input.f107=f107
    Input.ap=ap
    
    if alt > 500e3:
        gtd7d(Input, flags, output)
    else:
        gtd7(Input, flags, output)
    
    d = output.d
    t = output.t
    
    """
    DEFAULT OUTPUT VARIABLES:
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
    kb = 1.38064852e-23 # Boltzmann constant (m**2 kg)/(s**2 K)
    Na = 6.022140857e26 # avogadro number (molecules per kilomole)
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
    partial_p[0] = d[0]*kb*t[1] #He
    partial_p[1] = d[1]*kb*t[1] #O
    partial_p[2] = d[2]*kb*t[1] #N2
    partial_p[3] = d[3]*kb*t[1] #O2
    partial_p[4] = d[4]*kb*t[1] #AR
    partial_p[5] = d[6]*kb*t[1] #H
    partial_p[6] = d[7]*kb*t[1] #N
    partial_p[7] = d[8]*kb*t[1] #anomalous O

    #Assuming perfect gas, calculate atmospheric pressure
    pressure_mixture = np.sum(partial_p)
    
    temperature = t[1]    

    mole_fraction = np.divide(partial_p,pressure_mixture)
    
    molecular_weight_mixture = np.sum(np.multiply(mole_fraction,molecular_weights)) #kg/kmol
    
    mass_fractions = np.multiply(mole_fraction,
                                 np.divide(molecular_weights,molecular_weight_mixture))
    
    specific_gas_constants = R0/molecular_weights
    
    R_mixture = np.sum(np.multiply(mass_fractions,specific_gas_constants))
    
    number_density_mixture = np.sum(d) - d[5] 
    
    mean_free_path = (np.sqrt(2)*np.pi*4.15e-10**2*number_density_mixture)**-1
    
    eta = np.float64(1.458e-6*temperature**1.5/(temperature + 110.4)) # dynamic viscosity via sutherland law
    
    SoS = np.float64(np.sqrt(1.4*R_mixture*temperature))
    
    rho = d[5]

    return rho, pressure_mixture, temperature, R_mixture, mean_free_path, eta, molecular_weight_mixture, SoS
    
# US mutant Atmosphere
def US62_76(r,RE):
    """
    US62_76 is a very simple atmosphere model that uses the US76 standard 
    atmosphere below 80 km and the US62 standard atmosphere above 80km
    
    Inputs: 
        r: altitude 
        RE: radius of the Earth
    
    Outputs:
        rho: density 
        P: pressure
        T: temperature
        mfp: mean free path
        eta: viscosity (sutherland's law)
        MolW: molecular weight
        SoS: speed of sound  
    """
    #Some constants:
    #RE = 6378.137e3
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
    Z = (r - RE)*1e-3 # convert radius in m to altitude in km
    H = me2po(RE,Z) # geopotential altitude

    BLH = np.array([0., 11., 20., 32., 47., 51., 71., me2po(RE,86.), 
                    me2po(RE,100.), me2po(RE,110.), me2po(RE,120.), me2po(RE,150.), 
                    me2po(RE,160.), me2po(RE,170.), me2po(RE,190.), me2po(RE,230.), 
                    me2po(RE,300.), me2po(RE,400.), me2po(RE,500.), me2po(RE,600.), 
                    me2po(RE,700.)])

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

def me2po(RE,Z):
    """
    me2po converts geometric altitude to geopotential altitude -- the US
    standard atmosphere works in geopotential altitudes, which approximates the
    altitude of a pressure surface above the mean sea level.
    
    The reasoning for this is as follows: A change in geometric altitude will
    create a change in gravitational potential energy per unit mass (as the 
    effects of gravity become smaller as two objects move away from each other)
    
    Inputs:
        RE: Earth radius
        Z: Geometric altitude
    Outputs:
        H: Geopotential altitude
    """
    H = RE*Z/(RE + Z)
    return H