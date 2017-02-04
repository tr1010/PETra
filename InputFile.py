#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run Script: Script used to set inputs and run 6DOF Trajectory model

"""

def get_inputs():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #                            I N P U T S                            #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # Atmosphere options
    Atmosphere = 0 # 0 for default (currently my mutant US62/US76 atmosphere)
    
    # Aerodynamics solver options
    ContinuumMethod = 0 # 1 for Classic Newtonian, 0 for Modified Newtonian 
    
    # Spacecraft dimensions and parameters
    GeometryChoice = 0 # choice of box or plate or sphere (at the moment only box is available)
    scWidth = 1. # Width [m]
    scDepth = 1. # Depth [m]
    scHeight = 1. # Height [m]
    scDensity = 500. # Density of spacecraft [kg/m3]
  
    # Initial Spacecraft Conditions    
    FlightPathAngle = 2.21 # Positive down [deg]
    Heading = 102. # Measured clockwise from North [deg]
    Speed = 7875.8 # Velocity
    Latitude = 48.68  # [deg] Actually Declination. Geodetic Latitude input will be implemented soon
    Longitude = 123.68 # [deg]
    Altitude = 120.0e3 # Initial altitude [m]
    
    Pitch = 0. # Euler Angle 1 measured from Geocentric coordinate system
    Yaw = 0. # Euler Angle 2 measured from Geocentric coordinate system
    Roll = 0. # Euler Angle 3 measured from Geocentric coordinate system
    
    omega_x = 0. # Angular velocity about the BODY x- axis [deg/s]
    omega_y = 0. # Angular velocity about the BODY y- axis [deg/s]
    omega_z = 0. # Angular velocity about the BODY z- axis [deg/s]
    
    # Integration parameters
    ndt = 1001 # Number of output steps
    tmax = 500 # time [s]
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #          P O S T P R O C E S S I N G      O P T I O N S           #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # Postprocessing Options
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #               M O D E L      P A R A M E T E R S                  #
    #     (Don't mess with these unless you know what you're doing)     #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Earthmu = 398600.4415e9     # Earth gravity constant
    EarthRad = 6378.137e3       # Earth radius
    EarthJ2 = 1.0826230e-3      # Earth J2
    EarthOmega = 4.1780741e-3   # Earth angular velocity [deg/s]
    
    KnFM = 10       # Knudsen number for free-molecular flow
    KnCont  = 0.001 # Knudsen number for continuum flow
    a1 = 3./8       # Bridging function coefficient 1
    a2 = 1./8       # Bridging function coefficient 2
    
    SigN = 0.9      # Normal momentum accomodation coefficient
    SigT = 0.9      # Tangential momentum accomodation coefficient
    
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #               A S S E M B L E      I N P U T S                    #
    #   (Don't mess with this even if you do know what you're doing)    #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    inputDictionary = dict([('Atmosphere' , Atmosphere),
                       ('ContinuumMethod' , ContinuumMethod),
                       ('GeometryChoice' , GeometryChoice),
                       ('scWidth' , scWidth),
                       ('scDepth' , scDepth),
                       ('scHeight' , scHeight),
                       ('scDensity' , scDensity),
                       ('tmax' , tmax),
                       ('ndt' , ndt),
                       ('FlightPathAngle' , FlightPathAngle),
                       ('Heading' , Heading),
                       ('Speed' , Speed),
                       ('Latitude' , Latitude),
                       ('Longitude' , Longitude),
                       ('Altitude' , Altitude),
                       ('Pitch' , Pitch),
                       ('Yaw' , Yaw),
                       ('Roll' , Roll),
                       ('omega_x' , omega_x),
                       ('omega_y' , omega_y),
                       ('omega_z' , omega_z),
                       ('Earthmu' , Earthmu),
                       ('EarthRad' , EarthRad),
                       ('EarthJ2' , EarthJ2),
                       ('EarthOmega' , EarthOmega),
                       ('KnFM' , KnFM),
                       ('KnCont'  , KnCont),
                       ('a1' , a1),
                       ('a2' , a2),
                       ('SigN' , SigN),
                       ('SigT' , SigT)])

    return inputDictionary

















