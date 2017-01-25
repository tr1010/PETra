#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test Script for the Newtonian impact solver
Created on Mon Jan 23 14:22:01 2017

@author: tr1010
"""
import meshtools as mesh
import numpy as np
import AeroSolvers as Aero
import matplotlib.pyplot as plt

# Import geometry
CoG_off = np.zeros((3,1))
verts, surfs, normals, CoG = mesh.Cube(1.,CoG_off)

# Set free-stream conditions (sweep thru different attitudes?)
Ma = 20.
V_inf = np.array([Ma/6*np.sqrt(1.4*287.1*200.), Ma/5*np.sqrt(1.4*287.1*200.), Ma/2*np.sqrt(1.4*287.1*200.)])


# Calculate Pressure and force coefficients
Cp = Aero.NewtonSolver(normals, V_inf, Ma, 0)

# Plot results ?
