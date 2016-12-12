#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectories main
Created on Tue Nov 29 18:05:23 2016

@author: tr1010
"""
import numpy as np
import scipy as sp
import Trajectory as tr
import postprocessing as pp


def get_input():
    
    # spacecraft parameters
    #sc = np.array([2., 0.015, 2.2]) # m, A, Cd
    #sc = np.array([1140., 1.8*2.35, 2.2])
    sc = np.array([350., 4., 0.48])
    #sc = np.array([37.25, 0.5, 1.375])
    J = np.array([[80.,0,0],
                  [0, 120.0, 0],
                  [0, 0, 130.]])
    # time-stepping params:
    tmax = 10000
    ndt = 1001 
    t = np.linspace(0,tmax,ndt)
    
    
    # Earth constants based on WGS-84 Earth model: mu, RE, J2, ome
    earth = np.array([398600.4415e9, 6378.137e3, 1.0826230e-3, np.deg2rad(4.1780741e-3)])
    
    # Initial conditions: 
    #y0 = [120e3 + earth[1],np.deg2rad(19.73),np.deg2rad(-155.09), 7875., np.deg2rad(45), np.deg2rad(-2.2)]
    #y0 = [120e3 + earth[1],np.deg2rad(-79.8489),np.deg2rad(-10), 7835.7, np.deg2rad(100.), np.deg2rad(-0.051)] 
    y0 = np.array([6579.89967e3,np.deg2rad(-79.8489),np.deg2rad(-10), 7678.8, np.deg2rad(100.), np.deg2rad(0.55), 1, 0, 0, 0, 1., 0.1, -0.1]) # skipping?           
    #y0 = [120e3 + earth[1],np.deg2rad(48.68),np.deg2rad(123.68), 7875.8, np.deg2rad(282.), np.deg2rad(-2.21)]
    
    return earth, sc, t, y0, J

    

#Get input
earth, sc, t, y0, J = get_input()

# Integrate
sol = sp.integrate.odeint(tr.traj,y0,t,args=(sc, J, earth))

# Postprocess results
pp.plot(sol, sc, earth, t)