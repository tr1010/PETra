#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory model based on Equations of Motion in Chung & Lee
Created on Mon Nov 28 15:50:50 2016

@author: tr1010
"""
import numpy as np
import atmospheres as atmo

def traj(y, t, sc, earth):
    # Extract variables
    r, lat, lon, V, psi, gam= y
    mu, RE, J2, ome = earth
    m, A, Cd = sc
    if r <= earth[1]:
        dydt = np.zeros((6,))
    else:
        # Atmosphere calculation at new altitude
        rho, P, T, mfp, eta = atmo.US62_76(r)
    
        # Drag calculation
        D = 0.5*rho*V**2*A*Cd
        
        # Calculate trig functions
        slat = np.sin(lat)
        clat = np.cos(lat)
        tlat = np.tan(lat)
        sgam = np.sin(gam)
        cgam = np.cos(gam)
        spsi = np.sin(psi)
        cpsi = np.cos(psi)
    
        # Gravity equations
        gc = (mu/r**2)*(1 - 1.5*J2*(RE/r)**2*(3*slat**2 - 1))
        gd = -3*mu*J2*(RE/r**2)**2*slat*clat
        
        dydt = [V*sgam,
                V*cgam*cpsi/r,
                V*cgam*spsi/(r*clat),
                -D/m - gc*sgam + gd*cgam*cpsi - ome**2*r*clat*(cgam*cpsi*slat - sgam*clat),
                V*cgam*spsi*tlat/r - gd*spsi/(V*cgam) + ome**2*r*spsi*slat*clat/(V*cgam) - (2*ome/cgam)*(sgam*cpsi*clat - cgam*slat),
                V*cgam/r - gc*cgam/V - gd*sgam*cpsi/V + (ome**2*r*clat/V)*(sgam*cpsi*slat + cgam*clat) + 2*ome*spsi*clat]    
    

    return dydt