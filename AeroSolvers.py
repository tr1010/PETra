#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Newtonian impact theory panel method
Created on Mon Jan 23 10:18:06 2017

@author: tr1010
"""
import numpy as np
import math
    
def Aero_Calc(normals, Ma, Kn, RT, V, psi, gam, q1, q2, q3, q4, P, Q, R):
    """ 
    This function takes uses Newtonian impact theory and/or Schaaf & Chambre's 
    Molecular Flow model to calculate the Pressure coefficient on each panel of 
    the model geometry. It uses this pressure distributions to then calculate
    Aerodynamic Forces and Moments which are returned to be used in the 
    trajectory calculation
    """
    sgam = np.sin(gam)
    cgam = np.cos(gam)
    spsi = np.sin(psi)
    cpsi = np.cos(psi)
    
    # Assemble rotation matrices?
    Cw = np.array([[sgam, cgam*spsi, cgam*cpsi],
                   [0, cpsi, -spsi],
                   [-cgam, sgam*spsi, sgam*cpsi]])
    
    # In quaternion formulation?
    quat = np.array([q1,q2,q3,q4])
    
    S = np.array([[0., -q3, q2],
                  [q3, 0., -q1],
                  [-q2, q1, 0.]])
    
    temp = np.subtract(2*np.dot(quat,quat),np.multiply(2*q4,S))
    
    C = np.add((q4**2 - np.dot(quat,quat))*np.eye(3),temp)
    
    # Convert Quaternion to Euler angles -- not necessary?
#    Euler = np.array([np.arctan2(2*(q4*q1+q2*q3),1-2*(q1**2 + q2**2)),
#                      np.arcsin(2*(q4*q2 - q3*q1)),
#                      np.arctan2(2*(q4*q3 + q1*q2),1-2*(q2**2+q3**2))])
    
    # Find free stream velvec wrt to local horizon frame of reference
    temp = np.multiply(np.linalg.inv(Cw),np.linalg.inv(C))
    
    Vinf = np.multiply(temp,np.array([V,0.,0.]))
    
                      
    # Dependence on Knudsen number
    if Kn < 14.5:
        if Kn < 0.0146:
            # Calculate continuum pressure distribution
            Cp = NewtonSolver(normals,Vinf,Ma,0)
        else:
            # Calculate Transition pressure distribution
            #Cd = Cdcont + (Cdfm - Cdcont)*((1./3.)*np.log10(Kn/0.5) + 0.5113)
    else:
        # Calculate Free-molecular pressure distribution
        Cn, Ct = SchaafChambre(normals, Vinf, M, R, T, Tw, SigN, SigT)
        
    # Calculate aerodynamic forces by integrating over surface
    
    
    # Move back to local horizon frame of reference
    
    
    # Should treat inertia tensor here as well?
    
    
    
    #Calculate moments
    L = 0.
    M = 0.
    N = 0.
    # External torque around the centre of mass
    # This MUST be in the local horizon fixed frame of reference?
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
    return Cn, Ct
    
def CpMaxCalc(Ma):
    k = 1.4
    
    PO2_pinf = (((k+1)**2 * Ma**2)/(4*k*Ma**2 - 2*(k-1)))**(k/(1-k)) * \
               ((1-k+2*k*Ma**2)/(k+1))
               
    CpMax = (2/(k*Ma**2))*(PO2_pinf-1)
    
    return CpMax