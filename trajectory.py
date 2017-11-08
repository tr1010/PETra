#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory model based on Equations of Motion Karlgaard
Created on Mon Nov 28 15:50:50 2016

@author: tr1010 (Thomas Rees)

"""

import numpy as np
import atmospheres as atmo
import aerodynamics as aero
#from pprint import pprint
    
def traj_uvw(x, t, earth, mass, areas, normals, centroids, I, scLS, aero_params):
    """
    traj_uvw calculates the time derivative of the system state vector. See
    See documentation for a full description of the state vector and the 
    trajectory equations which are being solved, as well as the different frames
    of reference used.
    
    Inputs:
        x: state vector. Contains:
            r, phi, theta, u, v, w, e[0], e[1], e[2], e[3], angvel[0], 
            angvel[1], angvel[2] = x
        t: time variable (necessary for scipy.odeint)
        earth: python tuple of earth parameters: mu, RE, J2, omega
        mass: mass of spacecraft
        areas: n-element array of the spacecraft surface areas
        normals: 3xn element array of the outward-pointing unit normal vectors
        centroids: n-element array of the centroids of the spacecraft surfaces
        I: 3x3 spacecraft inertia tensor
        scLS: length-scale associated with spacecraft. By default, it is the 
              longest of the spacecraft's three dimensions
        aero_params: python tuple describing a number of parameters for the 
                     aerodynamics solver in the following order:
                     KnFM, KnCont, a1, a2, SigN, SigT = aero_params
            
    Outputs:
        dxdt: Time derivative of the state vector
    """
    # Unpack
    e = np.zeros(4)
    angvel = np.zeros(3)
    r, phi, theta, u, v, w, e[0], e[1], e[2], e[3], angvel[0], angvel[1], angvel[2] = x
    mu, RE, J2, ome = earth
    
    cphi = np.cos(phi)
    tphi = np.tan(phi)
    
    # Geodetic latitude to Declination transformation
    phigd = phi
    
    # Atmosphere calculation at new altitude and calculate useful aerodynamic
    # quantities
    #rho, P, T, mfp, eta, MolW, SoS = atmo.US62_76(r,earth[1])
    #R = 287.058
    rho, P, T, R, mfp, eta, MolW, SoS = atmo.nrlmsise00(172,0,29000,r-earth[1],phi,theta,16,150,150,4)
    
    Vinf = np.linalg.norm(np.array([u,v,w])) #speed of s/c
    Tw = 287.0
    Ma = Vinf/SoS
    Kn = mfp/scLS
    q_inf = 0.5*rho*Vinf**2
    
    if r-earth[1] <= 80e3 and Vinf/SoS < 3:
        # Check if Mach number falls below 5 (Newtonian theory fails) or
        # S/C hits ground
        # There is a better way of doing this which I should implement
        dxdt = np.zeros((13,)) 
    else:
     
        # Geocentric to body rotation matrix
        Gd = np.array([[e[0]**2 + e[1]**2 - e[2]**2 - e[3]**2, 2*(e[1]*e[2] + e[0]*e[3]), 2*(e[1]*e[3] - e[0]*e[2])],
                       [2*(e[1]*e[2] - e[0]*e[3]), e[0]**2 - e[1]**2 + e[2]**2 - e[3]**2, 2*(e[0]*e[1] + e[2]*e[3])],
                       [2*(e[1]*e[3] + e[0]*e[2]), 2*(e[2]*e[3] - e[0]*e[1]), e[0]**2 - e[1]**2 - e[2]**2 + e[3]**2]])
        
        #Geodetic to Geocentric transformation matrix
        Gphi = np.array([[np.cos(phi-phigd), 0, np.sin(phi-phigd)],
                         [0, 1, 0],
                         [-np.sin(phi-phigd), 0, np.cos(phi-phigd)]])
        
        #Geodetic to body rotation matrix
        Gcb = np.dot(Gd,Gphi)
    
        # freestream velocity in the body frame of reference
        VinfB = np.dot(Gcb,-np.array([u,v,w]))
    
        # Aerodynamics Calculations -- forces and moments in the body frame of reference
        AeroF, AeroM = aero.aero_calc(VinfB, areas, normals, centroids, Ma, Kn, R, T, q_inf, P, Tw, aero_params)
    
        # Transform Aerodynamic forces to Geocentric FoR
        AeroF_GC = np.dot(np.linalg.inv(Gcb),np.transpose(AeroF))

        # Solve eqns of rotational motion to Calculate angular velocity in body frame
        temp = np.transpose(np.subtract(AeroM,np.cross(angvel,np.dot(I,angvel))))
        angvel_dot = np.dot(np.linalg.inv(I),temp)
        
        # Unit Quaternion rotation eqns
        emat = np.array([[-e[1], -e[2], -e[3]],
                         [e[0], -e[3], e[2]],
                         [e[3], e[0], -e[1]],
                         [-e[2], e[1], e[0]]])
    
        edotrot = np.subtract(angvel,(1./r)*np.dot(Gcb,np.array([v, -u, -v*tphi])))
        edot = 0.5*np.dot(emat,edotrot)
        
        J = 1.5*J2*earth[1]**2
        # time derivative of state variables
        dxdt = [-w,
                u/r,
                v/(r*cphi) - ome,
                AeroF_GC[0]/mass + (u*w-v**2*tphi)/r - (mu*J/(r**4))*np.sin(2*phi),
                AeroF_GC[1]/mass + (u*v*tphi + v*w)/r,
                AeroF_GC[2]/mass - (u**2 + v**2)/r + mu/r**2 - (mu*J/(r**4))*(2-3*cphi**2),
                edot[0],
                edot[1],            
                edot[2],
                edot[3],
                angvel_dot[0],
                angvel_dot[1],
                angvel_dot[2]
                ]
            
    return dxdt
