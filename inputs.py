#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_inputs():
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #                            I N P U T S                            #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # Atmosphere options
    Atmosphere = 1 # 0 for mutant US62/US76 standard atmosphere, 1 for NRLMSISE-00
    
    #nrlmsise-00 options. 
    nrlmsise00_year = 0.
    nrlmsise00_doy = 172.
    nrlmsise00_sec = 29000.
    nrlmsise00_lst = 16.
    nrlmsise00_f107a = 150.
    nrlmsise00_f107 = 150.
    nrlmsise00_ap = 4.
    
    # Aerodynamics solver options
    ContinuumMethod = 0 # 1 for Classic Newtonian, 0 for Modified Newtonian 
    
    # Spacecraft dimensions and parameters
    GeometryChoice = 0 # choice of box or plate or sphere (at the moment only box is available)
    scWidth = 3.4 # Width [m]
    scDepth = 4.8 # Depth [m]
    scHeight = 2.35 # Height [m]
    scDensity = 30. # Density of spacecraft [kg/m3]
  
    # Initial Spacecraft Conditions    
    FlightPathAngle = 0. # Positive down [deg]
    Heading = 102. # Measured clockwise from North [deg]
    Speed = 7727. # Velocity
    Latitude = -79.8489  # [deg] Actually Declination. Geodetic Latitude input will be implemented soon
    Longitude = -10 # [deg]
    Altitude = 220e3 # Initial altitude [m]
    
    Pitch = 0. # Euler Angle 1 measured from Geocentric coordinate system
    Yaw = 0. # Euler Angle 2 measured from Geocentric coordinate system
    Roll = 0. # Euler Angle 3 measured from Geocentric coordinate system
    
    omega_x = 0. #3.14 # Angular velocity about the BODY x- axis [deg/s]
    omega_y = 0. #0.314 # Angular velocity about the BODY y- axis [deg/s]
    omega_z = 0. #314 # Angular velocity about the BODY z- axis [deg/s]
    
    # Integration parameters
    ndt = 1001 # Number of output steps
    tmax = 3000 # time [s]
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #          P O S T P R O C E S S I N G      O P T I O N S           #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # Postprocessing Options
    Run_name = 'Run01'
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #               M O D E L      P A R A M E T E R S                  #
    #                  (Not necessary to touch this)                    #
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
    #                 (Not necessary to touch this)                     #
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    inputDictionary = dict([('Atmosphere' , Atmosphere),
                            ('nrlmsise00_year' , nrlmsise00_year),
                            ('nrlmsise00_doy' , nrlmsise00_doy),
                            ('nrlmsise00_sec' , nrlmsise00_sec),
                            ('nrlmsise00_lst' , nrlmsise00_lst),
                            ('nrlmsise00_f107a' ,nrlmsise00_f107a),
                            ('nrlmsise00_f107' , nrlmsise00_f107)
                            ('nrlmsise00_ap' , nrlmsise00_ap)
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

















