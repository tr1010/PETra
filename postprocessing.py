#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:44:46 2016

@author: tr1010
"""
import numpy as np
import matplotlib.pyplot as plt
import atmospheres as atmo
import scipy as sp
import Trajectory as tr

def plot(sol, sc, earth, t):
    
    ndt = np.size(t)
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
    # temp dummy vars
    q1 = 0
    q2 = 0
    q3 = 0
    q4 = 0
    P = 0
    Q = 0
    R = 0
    for i in range(0,ndt):
        if sol[i,0] <= earth[1]:
            break
        else:   
            Cd[i], Ext = tr.Aero_Calc(Ma[i],Kn[i],V[i],RT[i], q1, q2, q3, q4, P, Q, R)
        
    
    decel = np.diff(V)/np.diff(t)/-9.81
    Qdot = 0.5*rho*V**3*sc[1]*Cd/20.
    soak = sp.integrate.cumtrapz(Qdot,t)
    frag = soak*0.5*rho[0:ndt-1]*V[0:ndt-1]**2
    KE = 0.5*sc[0]*V**2
    

    # Euler parameters to euler angles
#    pitch = 
#    roll = 
#    yaw = 
    # What about calculating deceleration?
    
    fig1 = plt.figure(1)
    plt.plot(t,h*1e-3)
    #fig1.suptitle('Altitude vs time')
    plt.xlabel('time [s]')
    plt.ylabel('altitude [km]')
    plt.savefig('Figures/altitude.png', bbox_inches='tight')
    plt.show
    
    fig2 = plt.figure(2)
    plt.plot(t,V*1e-3)
    #fig2.suptitle('Velocity vs time')
    plt.xlabel('time [s]')
    plt.ylabel('Velocity [km/s]')
    plt.savefig('Figures/Velocity.png', bbox_inches='tight')
    plt.show
    
    fig3 = plt.figure(3)
    plt.plot(long,lat)
    #fig3.suptitle('Latitude vs Longitude')
    plt.xlabel('Longitude [deg]')
    plt.ylabel('latitude [deg]')
    plt.savefig('Figures/LatLong.png', bbox_inches='tight')
    plt.show
    
    fig4 = plt.figure(4)
    plt.plot(t,Ma)
    #fig4.suptitle('Mach number vs time')
    plt.xlabel('time [s]')
    plt.ylabel('Mach Number [-]')
    plt.savefig('Figures/Mach.png', bbox_inches='tight')
    plt.show
    
    fig5 = plt.figure(5)
    plt.plot(t[0:ndt-1],decel)
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
    return 