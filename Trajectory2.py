#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory model based on Equations of Motion in Chung & Lee
Created on Mon Nov 28 15:50:50 2016

@author: tr1010
"""
import numpy as np
import atmospheres as atmo

def Cd_Calc(M,Kn,V,RT):
    # Dependence of continuum Cd on Mach number
    Mrange = np.array([0., 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 12., 100.])
    Cdcontrange = np.array([0.45, 0.5, 0.78, 1.0, 0.75, 0.62, 0.6, 0.55, 0.55])
    
    s = V/(2*RT)**0.5
    Cdfm = 1.75 + np.pi**0.5/(2*s)
    if Kn < 14.5:
        for i in range(0,9):
            if M > Mrange[i] and M <= Mrange[i+1]:
                Cdcont = Cdcontrange[i] + (Cdcontrange[i+1] - Cdcontrange[i])*(M - Mrange[i])/(Mrange[i+1] - Mrange[i])           
        if Kn < 0.0146:
            Cd = Cdcont
        else:
            Cd = Cdcont + (Cdfm - Cdcont)*((1./3.)*np.log10(Kn/0.5) + 0.5113)
    else:
        Cd = Cdfm
        
    return Cd

def traj(y, t, sc, earth):
    # Extract variables
    r, lat, lon, V, psi, gam = y
    mu, RE, J2, ome = earth
    m, A, Cd = sc
    if r <= earth[1]:
        dydt = np.zeros((6,))
    else:
        # Atmosphere calculation at new altitude
        rho, P, T, mfp, eta, MolW = atmo.US62_76(r)
        
        # Find Cd 
        Kn = mfp/A**0.5
        R = 8314.32/MolW
        RT = R*T
        SoS = (1.4*RT)**0.5
        Ma = V/SoS
        Cd = Cd_Calc(Ma,Kn,V,RT)        
        
        # Drag calculation
        D = 0.5*rho*V**2*A*Cd
        
        # Calculate trig functions
        slat = np.sin(lat)
        clat = np.cos(lat)
        tlat = np.tan(lat)
        sgam = np.sin(gam)
        cgam = np.cos(gam)
        tgam = np.tan(gam)
        spsi = np.sin(psi)
        cpsi = np.cos(psi)
    
        # Gravity equations
        phi = np.pi/2 - lat
        cphi = np.cos(phi)
        sphi = np.sin(phi)
        gc = mu*(1 - 1.5*J2*(3*cphi**2 - 1)*(RE/r)**2)/r**2
        gd = -3*mu*J2*sphi*cphi*(RE/r)*(RE/r)/r**2
        
        dydt = [V*sgam,
                V*cgam*cpsi/r,
                V*cgam*spsi/(r*clat),
                -D/m - gc*sgam + gd*cgam*cpsi - ome**2*r*clat*(cgam*cpsi*slat - sgam*clat),
                V*cgam*spsi*tlat/r - gd*spsi/(V*cgam) + ome**2*r*spsi*slat*clat/(V*cgam) - (2*ome)*(tgam*cpsi*clat - slat),
                (V/r - gc/V)*cgam - gd*sgam*cpsi/V + ome**2*r*clat*(sgam*cpsi*slat + cgam*clat)/V + 2*ome*spsi*clat]   
    

    return dydt