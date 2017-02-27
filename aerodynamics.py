#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aerodynamics.py includes the functions used to calculate the aerodynamic
moments and forces experienced by an object entering the atmosphere which are
required for uvw_trajectory function

Created on Mon Jan 23 10:18:06 2017

@author: tr1010 (Thomas Rees)
"""
import numpy as np
import math
    
def aero_calc(Vinf, areas, normals, centroids, Ma, Kn, R, T, q_inf, p_inf, Tw, aero_params):
    """ 
    aero_calc calculates the aerodynamic forces and moments acting on the
    spacecraft (in the spacecraft body frame of reference). In the case of 
    hypersonic continuum flow, aero_calc uses Newtonian impact theory to calculate 
    the pressure distribution over the body. In the free-molecular flow regime 
    aero_calc usesSchaaf & Chambre's formulation. For the transitional regime,
    the bridging function of Wilmoth et al. is used.
    
    Inputs:
        Vinf: a 3-vector describing the oncoming free-stream flow velocity in 
              the body frame of reference
        areas: a n-element array of the areas of each of the shape's surfaces.
        centroids: a 3xn numpy array of reals describing the position of each 
                   surface's centroid 
        Ma: Mach number of the free stream flow
        Kn: Knudsen number of the free stream flow   
        R: Gas constant of the free stream flow  
        T: Free-stream temperature of the flow
        q_inf: dynamic pressure of the free stream flow 
        p_inf: atmospheric pressure (static pressure of the free stream)
        Tw: Wall temperature (currently fixed to default 287 K)
        aero_params: python tuple describing a number of parameters for the 
                     aerodynamics solver in the following order:
                     KnFM, KnCont, a1, a2, SigN, SigT = aero_params
            
    
    Outputs:
        AeroF: a 3 vector describing the aerodynamic forces acting on the body
               in the body frame of reference
        AeroM: a 3 vector describing the aerodynamic moments acting on the body
               in the body frame of reference
    
    """
    
    numpans = np.size(normals,1)
    AeroF = np.zeros(3)
    AeroM = np.zeros(3)
    KnFM, KnCont, a1, a2, SigN, SigT = aero_params
    
    #calculate aerodynamic forces or moments
    # for continuum or transitional flow
    if Kn < KnFM:
        if Kn < KnCont:
            # Calculate continuum pressure distribution
            Cp = newton_solver(normals,Vinf,Ma,1)
            Ct = np.zeros(numpans)
            Pn = q_inf*Cp + p_inf
            St = q_inf*Ct
            # Sum across panels to calculate forces and moments
            for i in range(0,numpans):
                tempCont = -Pn[i]*normals[:,i] + St[i]*(np.cross(normals[:,i],np.cross(Vinf/np.linalg.norm(Vinf),normals[:,i])))  
                AeroF = AeroF + tempCont*areas[i]   
                AeroM = AeroM + (np.cross(centroids[:,i],tempCont))*areas[i]            
            
        else:
            # Calculate continuum & FM pressure dist if in transition
            CpCont = newton_solver(normals,Vinf,Ma,1)
            CtCont = np.zeros(numpans)
            CpFM, CtFM = schaaf_chambre(normals, Vinf, R, T, Tw, SigN, SigT)
            
            PnCont = q_inf*CpCont + p_inf
            StCont = q_inf*CtCont
            
            PnFM = q_inf*CpFM + p_inf
            StFM = q_inf*CtFM
            
            # Sum across panels to calculate forces and moments
            for i in range(0,numpans):
                tempCont = -PnCont[i]*normals[:,i] + StCont[i]*(np.cross(normals[:,i],np.cross(Vinf/np.linalg.norm(Vinf),normals[:,i])))  
                tempFM = -PnFM[i]*normals[:,i] + StFM[i]*(np.cross(normals[:,i],np.cross(Vinf/np.linalg.norm(Vinf),normals[:,i])))
                
                # as aero coefficients
                tempContCoeff = tempCont/(q_inf)
                tempFMCoeff = tempFM/(q_inf)
                
                # Apply bridging function (Wilmoth et al.)              
                ForceCoeff = tempContCoeff + (tempFMCoeff - tempContCoeff)*np.sin(np.pi*(a1 + a2*np.log10(Kn)))
                
                AeroF = AeroF + ForceCoeff*areas[i]*q_inf   
                AeroM = AeroM + (np.cross(centroids[:,i],ForceCoeff))*areas[i]

    # for free-molecular
    else:
        # Calculate Free-molecular pressure distribution
        CpFM, CtFM = schaaf_chambre(normals, Vinf, R, T, Tw, SigN, SigT)
        
        Pn = q_inf*CpFM + p_inf
        St = q_inf*CtFM
        
        for i in range(0,numpans):
                tempCont = -Pn[i]*normals[:,i] + St[i]*(np.cross(normals[:,i],np.cross(Vinf/np.linalg.norm(Vinf),normals[:,i])))  
                AeroF = AeroF + tempCont*areas[i]   
                AeroM = AeroM + (np.cross(centroids[:,i],tempCont))*areas[i]
                                          

    return AeroF, AeroM

def newton_solver(normals, Vinf, M, switch):
    """
    newton_solver calculates the pressure distribution on a body in a hypersonic
    flow using Newtonian impact theory.
    
    Inputs:
        normals:a 3xn vector of the outward pointing unit normal vectors for each
                of the surfaces making up the body.
        Vinf: a 3-vector describing the free-stream flow velocity in the body
              frame of reference
        M: free stream Mach number
        switch: boolean switch. 0 is Modified Newtonian, 1 is standard Newtonian
        
    Outputs:
        Cp: a n-element array of the pressure coefficients on each of the surfaces
            making up the body.
    """
    # declare/initialise variables
    totpans = np.size(normals,1)
    Cp = np.zeros((totpans,1))
    
    # Check if using modified Newtonian impact method or old-fashioned
    if switch == 1:
        CpMax = 2.
    else:
        CpMax = cp_max_calc(M)
    
    # Loop over panels to calculate pressure distribution
    for i in range(0,totpans):

        stheta = np.dot(Vinf,normals[:,i])/(np.linalg.norm(Vinf)*np.linalg.norm(normals[:,i]))
        
        # check if panel is shaded or not then calculate Cp
        if stheta < 0 and np.linalg.norm(stheta) > 1e-8:
            Cp[i] = CpMax*stheta**2.
        else:
            Cp[i] = 0
 
    return Cp
    
def schaaf_chambre(normals, Vinf, R, T, Tw, SigN, SigT):
    """
    schaaf_chambre uses Schaaf & Chambre's formulation to calculate the pressure
    distribution in a free-molecular flow regime.
    
    Inputs:
        normals:a 3xn vector of the outward pointing unit normal vectors for each
                of the surfaces making up the body.
        Vinf: a 3-vector describing the free-stream flow velocity in the body
              frame of reference
        R: Gas constant of the free stream flow  
        T: Free-stream temperature of the flow
        Tw: Wall temperature (currently fixed to default 287 K)
        SigN: Normal momentum accomodation coefficient
        SigT: tangential momentum accomodation coefficient
    
    Outputs:
        Cn: n-element array of the normal aerodynamic force coefficients
        Ct: n-element array of the tangential aerodynamic force coefficient
    """
    # declare/initialise variables
    totpans = np.size(normals,1)
    Cn = np.zeros((totpans,1))
    Ct = np.zeros((totpans,1))
    
    # Now calculate speed ratio and flow angle for each panel
    s = np.linalg.norm(Vinf)/(np.sqrt(2*R*T))
    
    # check if panel is shaded or not then calculate CN or Ct using Schaaf & Chambre
    for i in range(0,totpans):
        cdelta = np.dot(Vinf,normals[:,i])/(np.linalg.norm(Vinf)*np.linalg.norm(normals[:,i])) 
        sdelta = np.sqrt(1-cdelta**2)
        x = s*cdelta
        Gam1 = (x*np.exp(-x**2) + np.pi**0.5*(1+2*x**2)*(1+math.erf(x))/2)/np.pi**0.5
        Gam2 = (np.exp(-x**2) + np.pi**0.5*x*(1+math.erf(x)))/np.pi**0.5

        if cdelta < 0 and np.linalg.norm(cdelta) > 1e-8:
            Cn[i] = (1./s**2)*((2.-SigN)*Gam1 + SigN*((np.pi*Tw/T)**0.5)*Gam2/2.)
            Ct[i] = SigT*sdelta*Gam2/s
        else:
            Cn[i] = 0.
            Ct[i] = 0.

    return Cn, Ct
    
def cp_max_calc(Ma):
    """
    Calculates the maximum pressure coefficient for modified Newtonian flow
    
    Inputs:
        Ma: Free stream mach number
    Outputs:
        CpMax: Maximum pressure coefficient 
    """
    k = 1.4
    
    PO2_pinf = (((k+1)**2 * Ma**2)/(4*k*Ma**2 - 2*(k-1)))**(k/(1-k)) * \
               ((1-k+2*k*Ma**2)/(k+1))
               
    CpMax = (2/(k*Ma**2))*(PO2_pinf-1)
    
    return CpMax