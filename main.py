#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectories main
Created on Tue Nov 29 18:05:23 2016

@author: tr1010
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import atmospheres as atmo
import Trajectory2 as tr


# spacecraft parameters
sc = [350., 4., 2.2] # m, A, Cd
# Earth constants based on WGS-84 Earth model: mu, RE, J2, ome
earth = [398600.4415e9, 6378.137e3, 1.0826230e-3, np.deg2rad(4.1780741e-3)]

# Initial conditions: 
#y0 = [120e3 + earth[1],np.deg2rad(19.73),np.deg2rad(-155.09), 7875., np.deg2rad(45), np.deg2rad(-2.2)]
y0 = [6579.89967e3,np.deg2rad(-79.8489),np.deg2rad(-10), 7589.3, np.deg2rad(100.), np.deg2rad(0.55)]
#y0 = [6579.89967e3,np.deg2rad(-79.8489),np.deg2rad(-10), 7589.3, np.deg2rad(100.), np.deg2rad(0.55)]           
# Time array1000
t = np.linspace(0,5000,1001)

# Integrate
sol = sp.integrate.odeint(tr.traj,y0,t,args=(sc,earth))

# Postprocess results
h = sol[:,0] - earth[1]
V = sol[:,3]
gam = np.rad2deg(sol[:,5])

# Calculate Trajectory Flow characteristics
rho = np.zeros((1001,))
P = np.zeros((1001,))
T = np.zeros((1001,))
mfp = np.zeros((1001,))
eta = np.zeros((1001,))
Kn = np.zeros((1001,))
D = np.zeros((1001,))
SoS = np.zeros((1001,))
Ma = np.zeros((1001,))

for i in range(0,1001):
    if sol[i,0] <= earth[1]:
        break
    else:
        rho[i], P[i], T[i], mfp[i], eta[i] = atmo.US62_76(sol[i,0])
    
Kn = np.divide(mfp,sc[1]**0.5)
D = 0.5*rho*V**2*sc[1]*sc[2]
SoS = (1.4*287.*T)**0.5
Ma = np.divide(V,SoS)
plt.figure(1)
plt.plot(t,h)

plt.figure(2)
plt.plot(t,V)

plt.figure(4)
plt.plot(t,gam)


#Plot Results
#plt.figure(1)
#plt.plot(t,sol[:,0])
#plt.ylim((0,6500000))
#plt.show
#
#plt.figure(2)
#plt.plot(t,sol[:,5])
#plt.show


