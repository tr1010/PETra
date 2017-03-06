#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
function which processes the raw inputs
Created on Sat Feb  4 14:04:07 2017

@author: tr1010 (Thomas Rees)
"""
import numpy as np
import sys
import mesh_tools as mesh

def process_inputs(get_inputs):
    """
    process_inputs creates the individual variables necessary to run the main
    trajectory integrator. 
    
    Inputs:
        get_inputs: a python dictionary containing all of the inputs in the
                    inputs.py file. For the current version, the dictionary
                    contains 31 entries
    
    Outputs:
        earth: python tuple of earth parameters: mu, RE, J2, omega
        mass: mass of the spacecraft in kg
        areas: n element array containing the areas of each of the surfaces
               making up the spacecraft body
        normals: 3xn array containing the outward pointing unit normals of each
                 of the surfaces making up the spacecraft body
        centroids: 3xn array containing the locations of the centroids of each
                   of the surfaces making up the the spacecraft body with respect
                   to the spacecraft centre of mass
        I: 3x3 spacecraft inertia tensor
        t: array describing the times at which outputs are requested
        x0: the initial condition array for the scipy ode integrator Contains:
            r, phi, theta, u, v, w, e[0], e[1], e[2], e[3], angvel[0], angvel[1], 
            angvel[2] = x
        scLengthScale: length-scale associated with spacecraft. By default, it is the 
                       longest of the spacecraft's three dimensions
        aero_params: python tuple describing a number of parameters for the 
                     aerodynamics solver in the following order:
                     KnFM, KnCont, a1, a2, SigN, SigT = aero_params
        RawInputs: a python dictionary identical to get_inputs
    """
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
    earth = (RawInputs['Earthmu'], RawInputs['EarthRad'], 
                      RawInputs['EarthJ2'], RawInputs['EarthOmega']*d2r)  
    
    # Assemble initial condition vector
    t = np.linspace(0,RawInputs['tmax'],RawInputs['ndt'])
    
    x0 = np.array([RawInputs['Altitude'] + earth[1], RawInputs['Latitude']*d2r, 
                   RawInputs['Longitude']*d2r, u, v, w, e[0], e[1], e[2], e[3],
                   RawInputs['omega_x']*d2r, RawInputs['omega_y']*d2r, 
                   RawInputs['omega_z']*d2r])
    
    # Aerodynamics module parameters
    aero_params = (RawInputs['KnFM'], RawInputs['KnCont'], RawInputs['a1'],
                  RawInputs['a2'], RawInputs['SigN'], RawInputs['SigT'])

    return earth, mass, areas, normals, centroids, I, t, x0, scLengthScale, aero_params, RawInputs