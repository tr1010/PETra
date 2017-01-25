#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Newtonian impact theory panel method
Created on Mon Jan 23 10:18:06 2017

@author: tr1010
"""
import numpy as np
import math

def CoefficientCalculation(Cp, Areas, normals):
    totpans = np.size(Cp,1)
    forces = np.zeros((3,totpans))
    
   # Calculate local force vector on each panel
    for i in range(0,totpans):
        forces[:,i] = -Cp[i]*Areas[i]*normals[:,i]
    
    # Calculate Lift and Drag coefficients -- for this will need orientation/rotation matrices of s/c
    
    # Calculate Moment Coefficients -- need force coefficients along with orientation/rotation, CoG
    
    return forces
    
def Aero_Calc(Ma,Kn,V,RT, q1, q2, q3, q4, P, Q, R):
    # Dependence of continuum Cd on Mach number
    Mrange = np.array([0., 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 12., 100.])
    Cdcontrange = np.array([0.45, 0.5, 0.78, 1.0, 0.75, 0.62, 0.6, 0.55, 0.55])
    
    s = V/(2*RT)**0.5
    Cdfm = 1.75 + np.pi**0.5/(2*s)
    if Kn < 14.5:
        for i in range(0,9):
            if Ma > Mrange[i] and Ma <= Mrange[i+1]:
                Cdcont = Cdcontrange[i] + (Cdcontrange[i+1] - Cdcontrange[i])*(Ma - Mrange[i])/(Mrange[i+1] - Mrange[i])           
        if Kn < 0.0146:
            Cd = Cdcont
        else:
            Cd = Cdcont + (Cdfm - Cdcont)*((1./3.)*np.log10(Kn/0.5) + 0.5113)
    else:
        Cd = Cdfm
    
    L = 0.
    M = 0.
    N = 0.
    Ext = np.array([[L],
                   [M],
                   [N]])
    return Cd, Ext

def NewtonSolver(normals, Vinf, M, switch):   
    # declare/initialise variables
    totpans = np.size(normals,1)
    Cp = np.zeros((totpans,1))
    
    # Check if using modified Newtonian impact method or old-fashioned
    if switch == 1:
        CpMax = 2.
    else:
        CpMax = CpMaxCalc(M)
    
    # Loop over panels to calculate pressure distribution
    for i in range(0,totpans):

        stheta = np.dot(np.divide(Vinf,np.linalg.norm(Vinf)),normals[:,i])
        
        # check if panel is shaded or not then calculate Cp
        if stheta > 0:
            Cp[i] = CpMax*stheta**2.
        elif stheta < 0:
            Cp[i] = 0
        else:
            magnorm = np.linalg.norm(normals[:,i])
            magVinf = np.linalg.norm(Vinf)
            
            if (magnorm + magVinf) < np.max(np.array([magnorm,magVinf])):
                Cp[i] = CpMax*stheta**2.
            else:
                Cp[i] = 0.
 
    return Cp
    
def SchaafChambre(normals, Vinf, M, R, T, Tw, SigN, SigT):
    # declare/initialise variables
    totpans = np.size(normals,1)
    Cn = np.zeros((totpans,1))
    Ct = np.zeros((totpans,1))
    
    # Now calculate speed ratio and flow angle for each panel
    s = Vinf/(np.sqrt(2*R*T))
    
    # check if panel is shaded or not then calculate CN or Ct using Schaaf & Chambre
    for i in range(0,totpans):
        stheta = np.dot(np.divide(Vinf,np.linalg.norm(Vinf)),normals[:,i]) 
        ctheta = np.sqrt(1-stheta**2)
        if stheta > 0:
            Cn[i] = (1./s**2)*((((2-SigN)*s*stheta/np.sqrt(np.pi) + SigN*np.sqrt(Tw/T)/2.)*np.exp(-(s*stheta)**2)) + \
                    ((2-SigN)*((s*stheta)**2 + 0.5) + SigN*s*stheta*np.sqrt(np.pi * Tw/T)/2.) * \
                    (1+math.erf(s*stheta)))
            Ct[i] = (-SigT*ctheta/s/np.sqrt(np.pi))*(np.exp(-(s*stheta)**2) + np.sqrt(np.pi)*s*stheta*(1.+math.erf(s*stheta)))
        elif stheta < 0:
            Cn[i] = 0
            Ct[i] = 0
        else:
            magnorm = np.linalg.norm(normals[:,i])
            magVinf = np.linalg.norm(Vinf)
            
            if (magnorm + magVinf) < np.max(np.array([magnorm,magVinf])):
                Cn[i] = (1./s**2)*((((2-SigN)*s*stheta/np.sqrt(np.pi) + SigN*np.sqrt(Tw/T)/2.)*np.exp(-(s*stheta)**2)) + \
                    ((2-SigN)*((s*stheta)**2 + 0.5) + SigN*s*stheta*np.sqrt(np.pi * Tw/T)/2.) * \
                    (1+math.erf(s*stheta)))
                Ct[i] = (-SigT*ctheta/s/np.sqrt(np.pi))*(np.exp(-(s*stheta)**2) + np.sqrt(np.pi)*s*stheta*(1.+math.erf(s*stheta)))
            else:
                Cn[i] = 0
                Ct[i] = 0
    return
    
def CpMaxCalc(Ma):
    k = 1.4
    
    PO2_pinf = (((k+1)**2 * Ma**2)/(4*k*Ma**2 - 2*(k-1)))**(k/(1-k)) * \
               ((1-k+2*k*Ma**2)/(k+1))
               
    CpMax = (2/(k*Ma**2))*(PO2_pinf-1)
    
    return CpMax