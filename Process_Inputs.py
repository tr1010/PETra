#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 14:04:07 2017

@author: tr1010
"""
import numpy as np
import sys
import meshtools as mesh

def process_inputs(get_inputs):
    d2r = np.pi/180.
    
    RawInputs = get_inputs()
    
    # Do S/C geometry    
    BLengths = np.array([RawInputs['scWidth'],RawInputs['scHeight'],RawInputs['scDepth']])
    CoG_Off = np.array([0.,0.,0.]) # CoG Offset -- ignore for now
    scLengthScale = np.max(BLengths)
    
    if RawInputs['GeometryChoice'] == 0:
        verts, surfs, areas, normals, centroids, CoG, I, mass = mesh.Box(BLengths,RawInputs['scDensity'],CoG_Off)
    else:
        print('Error: Wrong Geometry Choice!')
        sys.exit()
        
    
    # Do S/C initial conditions    
    u = RawInputs['Speed']*np.cos(RawInputs['Heading']*d2r)*np.cos(RawInputs['FlightPathAngle']*d2r)
    v = RawInputs['Speed']*np.cos(RawInputs['FlightPathAngle']*d2r)*np.sin(RawInputs['Heading']*d2r)
    w = RawInputs['Speed']*np.sin(RawInputs['FlightPathAngle']*d2r)
    
    # Euler to quaternion
    sPitch_2 = np.sin(RawInputs['Pitch']*d2r/2.)
    cPitch_2 = np.cos(RawInputs['Pitch']*d2r/2.)
    sYaw_2 = np.sin(RawInputs['Yaw']*d2r/2.)
    cYaw_2 = np.cos(RawInputs['Yaw']*d2r/2.)
    sRoll_2 = np.sin(RawInputs['Roll']*d2r/2.)
    cRoll_2 = np.cos(RawInputs['Roll']*d2r/2.)

    e = np.zeros(4)
    e[0] = cRoll_2*cPitch_2*cYaw_2 + sRoll_2*sPitch_2*sYaw_2
    e[1] = sRoll_2*cPitch_2*cYaw_2 - cRoll_2*sPitch_2*sYaw_2
    e[2] = cRoll_2*sPitch_2*cYaw_2 + sRoll_2*cPitch_2*sYaw_2
    e[3] = cRoll_2*cPitch_2*sYaw_2 - sRoll_2*sPitch_2*cYaw_2
    
    
    # Earth constants based on WGS-84 Earth model: mu, RE, J2, ome
    earth = np.array([RawInputs['Earthmu'], RawInputs['EarthRad'], 
                      RawInputs['EarthJ2'], RawInputs['EarthOmega']*d2r])  
    
    # Assemble initial condition vector
    t = np.linspace(0,RawInputs['tmax'],RawInputs['ndt'])
    
    x0 = np.array([RawInputs['Altitude'] + earth[1], RawInputs['Latitude']*d2r, 
                   RawInputs['Longitude']*d2r, u, v, w, e[0], e[1], e[2], e[3],
                   RawInputs['omega_x']*d2r, RawInputs['omega_y']*d2r, 
                   RawInputs['omega_z']*d2r])

    return earth, mass, areas, normals, centroids, I, t, x0, scLengthScale, RawInputs