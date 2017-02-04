#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:44:46 2016

@author: tr1010
"""
import numpy as np
import matplotlib.pyplot as plt
import atmospheres as atmo

def postprocess(t, sol, earth, mass, areas, normals, centroids, I, x0, scLengthScale):
    r2d = 180./np.pi
    ndt = np.size(t)
    # Unpack solution array:
        
    #Position
    Altitude = (sol[:,0] - earth[1])*1e-3
    Latitude =  sol[:,1]*r2d # Actually Declination
    Longitude = sol[:,2]*r2d
    
    # Geocentric Velocity components
    u = sol[:,3]
    v = sol[:,4]
    w = sol[:,5]

    # Velocity components to Spherical Coords
    Speed = np.zeros(ndt)
    for i in range(0,ndt): 
        Speed[i] = np.linalg.norm(np.array([u[i],v[i],w[i]]))
        
    FlightPathAngle = np.arcsin(np.divide(w,Speed))*r2d
    Heading = np.arccos(np.divide(u,Speed*np.cos(FlightPathAngle/r2d)))*r2d
    
    # Quaternion to Euler Angles
    Pitch = np.rad2deg(-np.arcsin(2*(sol[:,7]*sol[:,9] - sol[:,8]*sol[:,6])))
    Yaw = np.rad2deg(np.arctan2(2*(sol[:,8]*sol[:,7]+sol[:,6]*sol[:,9]),sol[:,6]**2+sol[:,7]**2 - sol[:,8]**2-sol[:,9]**2))
    Roll = np.rad2deg(np.arctan2(2*(sol[:,8]*sol[:,9]+sol[:,6]*sol[:,7]),sol[:,6]**2-sol[:,7]**2 - sol[:,8]**2+sol[:,9]**2))
    
    # Angular Velocities (around body axes)
    wx = sol[:,10]
    wy = sol[:,11]
    wz = sol[:,12]

    # Now calculate interesting outputs -- 
    # Mach number, Knudsen number, Energy, Force, Dynamic Pressure
    
    # Get atmospheric quantities at each point (this should be in main function
    # I promise I'll fix it soon)
    rho = np.zeros(ndt)
    P = np.zeros(ndt)
    T = np.zeros(ndt)
    mfp = np.zeros(ndt)
    eta = np.zeros(ndt)
    MolW = np.zeros(ndt)
    SoS = np.zeros(ndt)
    
    for i in range(0,ndt):
        if sol[i,0] <= earth[1]:
            break
        else:
            rho[i], P[i], T[i], mfp[i], eta[i], MolW[i], SoS[i]  = atmo.US62_76(sol[i,0])
    
    Mach = np.divide(Speed,SoS)
    Kn = scLengthScale/mfp
    Force = np.diff(Speed)/np.diff(t)/-9.81
    DynamicPressure = 0.5*np.multiply(rho,np.square(Speed))
    Energy = 0.5*mass*np.square(Speed) + mass*9.81*Altitude*1e3
    
    # Put all results in dictionary
    pp_output = dict([('FlightPathAngle' , FlightPathAngle),
                      ('Heading' , Heading),
                      ('Speed' , Speed),
                      ('Latitude' , Latitude),
                      ('Longitude' , Longitude),
                      ('Altitude' , Altitude),
                      ('Pitch' , Pitch),
                      ('Yaw' , Yaw),
                      ('Roll' , Roll),
                      ('omega_x' , wx),
                      ('omega_y' , wy),
                      ('omega_z' , wz),
                      ('Mach' , Mach),
                      ('Knudsen' , Kn),
                      ('Force' , Force),
                      ('DynamicPressure' , DynamicPressure),
                      ('Energy' , Energy)])
    
    # Plot Results in Python for user sanity check
    plotfn(0,t,pp_output['FlightPathAngle'], 'time', 'Flight Path Angle', 's', 'deg')
    plotfn(1,t,pp_output['Heading'], 'time', 'Heading', 's', 'deg')
    plotfn(2,t,pp_output['Speed'], 'time', 'Speed', 's', 'm/s')
    plotfn(3,pp_output['Longitude'],pp_output['Latitude'], 'Longitude', 'Latitude', 'deg', 'deg')
    plotfn(4,t,pp_output['Altitude'], 'time', 'Altitude', 's', 'km')
    plotfn(5,t,pp_output['Mach'], 'time', 'Mach Number', 's', '-')
    plotfn(6,t,pp_output['Knudsen'], 'time', 'Knudsen number', 's', '-')
    plotfn(7,t[0:ndt-1],pp_output['Force'], 'time', 'Force', 's', 'g')
    plotfn(8,t,pp_output['DynamicPressure']*1e-3, 'time', 'Dynamic Pressure', 's', 'kPa')
    plotfn(9,t,pp_output['Energy']*1e-9, 'time', 'Energy', 's', 'GJ')
    
    plt.figure(10)
    fig,ax = plt.subplots()
    plt.title('Euler angles')
    ax.plot(t,pp_output['Pitch'],label='pitch')
    ax.plot(t,pp_output['Yaw'], label='yaw')
    ax.plot(t,pp_output['Roll'],label = 'roll')
    legend = ax.legend()
    plt.xlabel('time [s]')
    plt.ylabel('angle [deg]')
    plt.show()
    
    return pp_output
    
    
    
def plotfn(figno, x, y, xname, yname, xunits, yunits):
    
    fig = plt.figure(figno)
    plt.plot(x,y)
    fig.suptitle(xname + ' vs ' + yname)
    plt.xlabel(xname + ' [' + xunits +']')
    plt.ylabel(yname + ' [' + yunits +']')
    plt.savefig('Figures/'+ xname + '-' + yname+'.png', bbox_inches='tight')
    plt.show
    
    return