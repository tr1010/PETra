#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Newtonian impact theory panel method
Created on Mon Jan 23 10:18:06 2017

@author: tr1010
"""
import numpy as np
import math
    
def Aero_Calc(Vinf, areas, normals, centroids, Ma, Kn, R, T, q_inf, p_inf, Tw):
    """ 
    This function takes uses Newtonian impact theory and/or Schaaf & Chambre's 
    Molecular Flow model to calculate the Pressure coefficient on each panel of 
    the model geometry. It uses this pressure distributions to then calculate
    Aerodynamic Forces and Moments which are returned to be used in the 
    trajectory calculation (in the body frame of reference)
    """
    
    numpans = np.size(normals,1)
    AeroF = np.zeros(3)
    AeroM = np.zeros(3)
    
    #calculate aerodynamic forces or moments
    # for continuum or transitional flow
#    if Kn < 10:
#        if Kn < 0.001:
#            # Calculate continuum pressure distribution
#            Cp = NewtonSolver(normals,Vinf,Ma,0)
#            Ct = np.zeros(numpans)
#            Pn = q_inf*Cp + p_inf
#            St = q_inf*Ct
#            # Sum across panels to calculate forces and moments
#            for i in range(0,numpans):
#                tempCont = -Pn[i]*normals[:,i] + St[i]*(np.cross(normals[:,i],np.cross(Vinf,normals[:,i])))  
#                AeroF = AeroF + tempCont*areas[i]   
#                AeroM = AeroM + (np.cross(centroids[:,i],tempCont))*areas[i]            
#            
#        else:
#            # Calculate continuum & FM pressure dist if in transition
#            CpCont = NewtonSolver(normals,Vinf,Ma,0)
#            CtCont = np.zeros(numpans)
#            CpFM, CtFM = SchaafChambre(normals, Vinf, Ma, R, T, Tw, SigN = 0.92, SigT = 0.92)
#            
#            PnCont = q_inf*CpCont + p_inf
#            StCont = q_inf*CtCont
#            
#            PnFM = q_inf*CpFM + p_inf
#            StFM = q_inf*CtFM
#            
#            # Sum across panels to calculate forces and moments
#            for i in range(0,numpans):
#                tempCont = -PnCont[i]*normals[:,i] + StCont[i]*(np.cross(normals[:,i],np.cross(Vinf,normals[:,i])))  
#                tempFM = -PnFM[i]*normals[:,i] + StFM[i]*(np.cross(normals[:,i],np.cross(Vinf,normals[:,i])))
#                
#                # as aero coefficients
#                tempContCoeff = tempCont/(q_inf)
#                tempFMCoeff = tempFM/(q_inf)
#                
#                # Apply bridging function (Wilmoth et al.)
#                a1 = 3./8
#                a2 = 1./8                
#                ForceCoeff = tempContCoeff + (tempFMCoeff - tempContCoeff)*np.sin(np.pi*(a1 + a2*np.log10(Kn)))
#                
#                AeroF = AeroF + ForceCoeff*areas[i]*q_inf   
#                AeroM = AeroM + (np.cross(centroids[:,i],ForceCoeff))*areas[i]
#                                          
#                #Cd = Cdcont + (Cdfm - Cdcont)*((1./3.)*np.log10(Kn/0.5) + 0.5113)
#    
#    # for free-molecular
#    else:
#        # Calculate Free-molecular pressure distribution
#        CpFM, CtFM = SchaafChambre(normals, Vinf, Ma, R, T, Tw, SigN = 0.92, SigT = 0.92)
#        
#        Pn = q_inf*CpFM + p_inf
#        St = q_inf*CtFM
#        
#        for i in range(0,numpans):
#                tempCont = -Pn[i]*normals[:,i] + St[i]*(np.cross(normals[:,i],np.cross(Vinf,normals[:,i])))  
#                AeroF = AeroF + tempCont*areas[i]   
#                AeroM = AeroM + (np.cross(centroids[:,i],tempCont))*areas[i]
                                          
    
    Cp = NewtonSolver(normals,Vinf,Ma,1)
    Ct = np.zeros(numpans)
    Pn = q_inf*Cp + p_inf
    St = q_inf*Ct
    # Sum across panels to calculate forces and moments
    for i in range(0,numpans):
        tempCont = -Pn[i]*normals[:,i] + St[i]*(np.cross(normals[:,i],np.cross(Vinf,normals[:,i])))  
        AeroF = AeroF + tempCont*areas[i]   
        AeroM = AeroM + (np.cross(centroids[:,i],tempCont))*areas[i] 

    return AeroF, AeroM

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

        stheta = np.dot(Vinf,normals[:,i])/(np.linalg.norm(Vinf)*np.linalg.norm(normals[:,i]))
        
        # check if panel is shaded or not then calculate Cp
        if stheta < 0:
            Cp[i] = CpMax*stheta**2.
        elif stheta > 0:
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
    s = np.linalg.norm(Vinf)/(np.sqrt(2*R*T))
    
    # check if panel is shaded or not then calculate CN or Ct using Schaaf & Chambre
    for i in range(0,totpans):
        stheta = -np.dot(np.divide(Vinf,np.linalg.norm(Vinf)),normals[:,i]) 
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