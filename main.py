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
import atmospheres as atmo

def get_input():
    # geom =[m, areas, normals, centroids, J]
    BLengths = np.array([3.4,1.8,2.35])
    CoG_Off = np.array([0.2,0.,0.])
    verts, surfs, areas, normals, centroids, CoG, I, mass = mesh.Box(BLengths,800,CoG_Off)
    scLS = np.max(BLengths)
    #geom = (mass, areas, normals, centroids, I)
    # time-stepping params:
    tmax = 300
    ndt = 1001 
    t = np.linspace(0,tmax,ndt)
    
    # Earth constants based on WGS-84 Earth model: mu, RE, J2, ome
    earth = np.array([398600.4415e9, 6378.137e3, 1.0826230e-3, np.deg2rad(4.1780741e-3)])   
    
    # Initial conditions: 
    # x0 = [r, theta, phi, u, v, w, e0, e1, e2, e3, angvel1, angvel2, angvel3]
    # Skipping trajectory
#    FPA = 0.54681217*np.pi/180.
#    Head = 99.955734*np.pi/180.
#    vel = 7.58930433867e3 # if this is too low calculation of CpMax fails
#    lat = -79.8489182889*np.pi/180.
#    long = -10*np.pi/180.0
#    alt = 405.0e3
    
    FPA = 2.21*np.pi/180.
    Head = (282.-180)*np.pi/180.
    vel = 7875.8 # if this is too low calculation of CpMax fails
    lat = 48.68*np.pi/180.
    long = 123.68*np.pi/180.
    alt = 120.0e3
    
    u = vel*np.cos(Head)*np.cos(FPA)
    v = vel*np.cos(FPA)*np.sin(Head)
    w = vel*np.sin(FPA)
    
    # initial orientation and angular velocities
    e = np.zeros(4)
    e[0] = 1.
    angvel = np.zeros(3)#np.array([0.005,0.005,0.005])
    
    x0 = np.array([alt + earth[1], lat, long, u, v, w, e[0], e[1], e[2], e[3], angvel[0], angvel[1], angvel[2]])
    
    return earth, mass, areas, normals, centroids, I, t, x0, scLS

    

#Get input
earth, mass, areas, normals, centroids, I, t, x0, scLS = get_input()

# test
#dxdt = tr.traj_uvw(x0, t, earth, mass, areas, normals, centroids, I)
# Integrate
sol = sp.integrate.odeint(tr.traj_uvw,x0,t,args=(earth, mass, areas, normals, centroids, I, scLS))




# postprocessing

ndt = np.size(t)
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
Vel = np.zeros(ndt,)       
for i in range(0,ndt):
    if sol[i,0] <= earth[1]:
        break
    else:
            rho[i], P[i], T[i], mfp[i], eta[i], MolW[i], SoS[i] = atmo.US62_76(sol[i,0])
            
    Vel[i] = np.linalg.norm(np.array([sol[i,3],sol[i,4],sol[i,5]]))

Ma = np.divide(Vel,SoS)


# quaternion to euler angles
pitch = np.rad2deg(-np.arcsin(2*(sol[:,7]*sol[:,9] - sol[:,8]*sol[:,6])))
yaw = np.rad2deg(np.arctan2(2*(sol[:,8]*sol[:,7]+sol[:,6]*sol[:,9]),sol[:,6]**2+sol[:,7]**2 - sol[:,8]**2-sol[:,9]**2))
roll = np.rad2deg(np.arctan2(2*(sol[:,8]*sol[:,9]+sol[:,6]*sol[:,7]),sol[:,6]**2-sol[:,7]**2 - sol[:,8]**2+sol[:,9]**2))

speed = np.zeros(np.size(sol,0))
for i in range(0,np.size(sol,0)):
    speed[i] = np.sqrt(sol[i,3]**2+sol[i,4]**2+sol[i,5]**2)
    
plt.figure()
plt.title('speed')
plt.plot(t,speed)
plt.figure()
plt.title('mach')
plt.plot(t,Ma)
plt.figure()
plt.title('altitude')
plt.plot(t,(sol[:,0]-earth[1])*1e-3)
fig,ax = plt.subplots()
plt.title('Euler angles')
ax.plot(t,pitch,label='pitch')
ax.plot(t,yaw, label='yaw')
ax.plot(t,roll,label = 'roll')
legend = ax.legend()
plt.show()
plt.figure()
plt.title('lat-long: something wrong?!??!')
plt.plot(np.rad2deg(sol[:,2]),np.rad2deg(sol[:,1]))

plt.figure()
plt.title('Mach-altitude')
plt.plot(Ma,(sol[:,0]-earth[1])*1e-3)


# Postprocess results
#pp.plot(sol, sc, earth, t)