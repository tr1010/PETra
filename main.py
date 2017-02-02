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
import meshtools as mesh
import matplotlib.pyplot as plt

def get_input():
    # geom =[m, areas, normals, centroids, J]
    verts, surfs, areas, normals, centroids, CoG, I, mass = mesh.Box(1,800,np.array([0.,0.,0.]))
    #geom = (mass, areas, normals, centroids, I)
    # time-stepping params:
    tmax = 1000
    ndt = 1001 
    t = np.linspace(0,tmax,ndt)
    
    # Earth constants based on WGS-84 Earth model: mu, RE, J2, ome
    earth = np.array([398600.4415e9, 6378.137e3, 1.0826230e-3, np.deg2rad(4.1780741e-3)])   
    
    # Initial conditions: 
    # x0 = [r, theta, phi, u, v, w, e0, e1, e2, e3, angvel1, angvel2, angvel3]
    FPA = 0.54681217*np.pi/180.
    Head = 99.955734*np.pi/180.
    vel = 7.58930433867e3
    lat = -79.8489182889*np.pi/180.
    long = -10*np.pi/180.0
    
#    FPA = 13*np.pi/180.
#    Head = 77.7*np.pi/180.
#    vel = 5500. #7.58930433867e2
#    lat = 69.36*np.pi/180.
#    long = 19.77*np.pi/180.0
    
    u = vel*np.cos(Head)*np.cos(FPA)
    v = vel*np.cos(FPA)*np.sin(Head)
    w = vel*np.sin(FPA)
    
    # initial orientation and angular velocities
    e = np.zeros(4)
    e[0] = 1.
    angvel = np.zeros(3)
    
    x0 = np.array([6579.89967e3, lat, long, u, v, w, e[0], e[1], e[2], e[3], angvel[0], angvel[1], angvel[2]])
    
    return earth, mass, areas, normals, centroids, I, t, x0

    

#Get input
earth, mass, areas, normals, centroids, I, t, x0 = get_input()

# test
#dxdt = tr.traj_uvw(x0, t, earth, mass, areas, normals, centroids, I)
# Integrate
sol = sp.integrate.odeint(tr.traj_uvw,x0,t,args=(earth, mass, areas, normals, centroids, I))

# quaternion to euler angles
pitch = np.rad2deg(-np.arcsin(2*(sol[:,7]*sol[:,9] - sol[:,8]*sol[:,6])))
yaw = np.rad2deg(np.arctan2(2*(sol[:,8]*sol[:,7]+sol[:,6]*sol[:,9]),sol[:,6]**2+sol[:,7]**2 - sol[:,8]**2-sol[:,9]**2))
roll = np.rad2deg(np.arctan2(2*(sol[:,8]*sol[:,9]+sol[:,6]*sol[:,7]),sol[:,6]**2-sol[:,7]**2 - sol[:,8]**2+sol[:,9]**2))

speed = np.zeros(np.size(sol,0))
for i in range(0,np.size(sol,0)):
    speed[i] = np.sqrt(sol[i,3]**2+sol[i,4]**2+sol[i,5]**2)
    
plt.figure()
plt.plot(t,speed)
plt.figure()
plt.plot(t,sol[:,0]-earth[1])
# Postprocess results
#pp.plot(sol, sc, earth, t)