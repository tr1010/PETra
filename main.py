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
#sc = np.array([2., 0.015, 2.2]) # m, A, Cd
#sc = np.array([1140., 1.8*2.35, 2.2])
sc = np.array([350., 4., 0.48])
#sc = np.array([37.25, 0.5, 1.375])

# time-stepping params:
tmax = 10000
ndt = 1001 
t = np.linspace(0,tmax,ndt)


# Earth constants based on WGS-84 Earth model: mu, RE, J2, ome
earth = np.array([398600.4415e9, 6378.137e3, 1.0826230e-3, np.deg2rad(4.1780741e-3)])

# Initial conditions: 
#y0 = [120e3 + earth[1],np.deg2rad(19.73),np.deg2rad(-155.09), 7875., np.deg2rad(45), np.deg2rad(-2.2)]
#y0 = [120e3 + earth[1],np.deg2rad(-79.8489),np.deg2rad(-10), 7835.7, np.deg2rad(100.), np.deg2rad(-0.051)] 
y0 = [6579.89967e3,np.deg2rad(-79.8489),np.deg2rad(-10), 7678.8, np.deg2rad(100.), np.deg2rad(0.55)] # skipping?           
#y0 = [120e3 + earth[1],np.deg2rad(48.68),np.deg2rad(123.68), 7875.8, np.deg2rad(282.), np.deg2rad(-2.21)]

# Integrate
sol = sp.integrate.odeint(tr.traj,y0,t,args=(sc,earth))

# Postprocess results
h = sol[:,0] - earth[1]
lat = np.rad2deg(sol[:,1])
long = np.rad2deg(sol[:,2])
V = sol[:,3]
psi = np.rad2deg(sol[:,4])
gam = np.rad2deg(sol[:,5])

# Calculate Trajectory Flow characteristics
rho = np.zeros((ndt,))
P = np.zeros((ndt,))
T = np.zeros((ndt,))
mfp = np.zeros((ndt,))
eta = np.zeros((ndt,))
Kn = np.zeros((ndt,))
D = np.zeros((ndt,))
SoS = np.zeros((ndt,))
Ma = np.zeros((ndt,))
MolW = np.zeros((ndt,))
Re = np.zeros((ndt,))
Cd = np.zeros(ndt,)
for i in range(0,ndt):
    if sol[i,0] <= earth[1]:
        break
    else:
        rho[i], P[i], T[i], mfp[i], eta[i], MolW[i] = atmo.US62_76(sol[i,0])

Kn = np.divide(mfp,sc[1]**0.5)
D = 0.5*rho*V**2*sc[1]*Cd
SoS = (1.4*287.*T)**0.5
Ma = np.divide(V,SoS)
Re = rho*V*sc[1]**0.5/eta
R = 8314.32/MolW
RT = R*T

for i in range(0,ndt):
    if sol[i,0] <= earth[1]:
        break
    else:   
        Cd[i] = tr.Cd_Calc(Ma[i],Kn[i],V[i],RT[i])
    

decel = np.diff(V)/np.diff(t)/-9.81
Qdot = 0.5*rho*V**3*sc[1]*Cd/20.
soak = sp.integrate.cumtrapz(Qdot,t)
frag = soak*0.5*rho[0:1000]*V[0:1000]**2
KE = 0.5*sc[0]*V**2

# What about calculating deceleration?

fig1 = plt.figure(1)
plt.plot(t,h*1e-3)
fig1.suptitle('Altitude vs time')
plt.xlabel('time [s]')
plt.ylabel('altitude [km]')
plt.show

fig2 = plt.figure(2)
plt.plot(t,V*1e-3)
fig2.suptitle('Velocity vs time')
plt.xlabel('time [s]')
plt.ylabel('Velocity [km/s]')
plt.show

fig3 = plt.figure(3)
plt.plot(long,lat)
fig3.suptitle('Latitude vs Longitude')
plt.xlabel('Longitude [deg]')
plt.ylabel('latitude [deg]')
plt.show

fig4 = plt.figure(4)
plt.plot(t,Ma)
fig4.suptitle('Mach number vs time')
plt.xlabel('time [s]')
plt.ylabel('Mach Number [-]')
plt.show

fig5 = plt.figure(5)
plt.plot(t[0:1000],decel)
fig5.suptitle('deceleration vs time')
plt.xlabel('time [s]')
plt.ylabel('deceleration [g]')
plt.show

fig6 = plt.figure(6)
plt.plot(t,Qdot*1e-6)
fig6.suptitle('Heat flux vs time')
plt.xlabel('time [s]')
plt.ylabel('Heat flux [MW]')
plt.show

fig7 = plt.figure(7)
plt.plot(frag*1e-12,h[0:1000]*1e-3)
fig7.suptitle('Fragmentation criterion vs height')
plt.ylabel('height [km]')
plt.xlabel('FC = Qdot * qinf [TW.Pa]')
plt.show

fig8 = plt.figure(8)
plt.plot(t,KE*1e-9,'b')
fig8.suptitle('Kinetic Energy vs time')
plt.xlabel('time [s]')
plt.ylabel('Kinetic Energy [GJ]')
plt.show

fig9 = plt.figure(9)
plt.plot(t,0.0005*rho*V**2)
fig9.suptitle('Dynamic Pressure vs time')
plt.xlabel('time [s]')
plt.ylabel('dynamic pressure [kPa]')
plt.show

#fig10 = plt.figure(10)
#plt.plot(t,Qdot*1e-6)
#fig10.suptitle('Heat flux vs time')
#plt.xlabel('time [s]')
#plt.ylabel('Heat flux [MW]')
#plt.show