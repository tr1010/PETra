#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectory model based on Equations of Motion in Karlgaard et al.
Created on Tue Nov 22 17:11:27 2016

@author: tr1010

a contains aerodynamic accelerations:
    ax, ay, az
earth contains WGS-84 Earth model constants:
    earth = Radius, eccentricity, mu, J2, omega
    
Note:
    phi is DECLINATION
    theta is LONGITUDE

"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def traj(y, t, a, earth):
    r, phi, theta, u, v, w = y
    mu = earth[2]
    J2 = earth[3]
    Omega = earth[4]
    dydt = [-w, 
            u/r, 
            v/(r*np.cos(phi)) - Omega, 
            a[0] + (1/r)*(u*w-v*v*np.tan(phi)) - (3*mu*J2/2*r**4)*np.sin(2*phi),
            a[1] + (1/r)*(u*v*np.tan(phi) + v*w),
            a[2] - (1/r)*(u*u + v*v) + mu/(r*r) - (3*mu*J2/2*r**4)*(2-3*np.cos(phi)**2)]
    
    return dydt
    
    
# constants
earth = [6378.137e3, 
         8.1819191e-2, 
         398600.4415e9, 
         1.0826230e-3, 
         np.deg2rad(4.1780741e-3)]
a = [-200, -200, 0]
#Velocity 7875, entry angle 0.1 deg, initial alt 120km, latlong as below
y0 = [120000e3 + earth[0], 
      np.deg2rad(19.73), 
      np.deg2rad(-155.09), 
      3941.12, 
      3941.12,
      -14]
      
t = np.linspace(0,5000,5001)

sol = sp.integrate.odeint(traj,y0,t,args=(a,earth))
plt.figure(1)
plt.plot(t,sol[:,0])
plt.ylim((0,6500000))
plt.show

plt.figure(2)
plt.plot(t,sol[:,5])
plt.show