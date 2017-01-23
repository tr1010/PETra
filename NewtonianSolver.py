#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Newtonian impact theory panel method
Created on Mon Jan 23 10:18:06 2017

@author: tr1010
"""
import numpy as np
import meshlibrary as ml

def CalcNormals(verts, surfs):
    
    totpans = np.size(surfs, 1)
    normals = np.zeros((3,totpans))
    
    for i in range(0,totpans):
        # Extract points making up the panel
        v1 = verts[:,surfs[1,i]] - verts[:,surfs[0,i]]
        v2 = verts[:,surfs[3,i]] - verts[:,surfs[0,i]]
        normals[:,i] = np.cross(v1, v2) # need to think about connectivity and the direction 'around' the surface
        
    return normals
    
def CpMaxCalc(Ma):
    k = 1.4
    
    PO2_pinf = (((k+1)**2 * Ma**2)/(4*k*Ma**2 - 2*(k-1)))**(k/(1-k)) 
    * ((1-k+2*k*Ma**2)/(k+1))
    
    CpMax = (2/(k*Ma**2))*(PO2_pinf-1)
    
    return CpMax

def NewtonSolver(normals, Vinf, M):
    
    totpans = np.size(normals,1)
    CpMax = CpMaxCalc(Ma)
    
    Cp = np.zeros((1,totpans))
    for i in range(0,totpans):

        stheta = np.dot(np.divide(Vinf,np.linalg.norm(Vinf)),normals[:,i])
        # check if panel is shaded or not
        #if normals[:,i] + np.divide(Vinf,np.linalg.norm(Vinf)) < normals[:,i]
        Cp[i] = CpMax*stheta**2
        
    return Cp
    
def CoefficientCalculation(Cp, Areas):
    totpans = np.size(Cp,1)
    forces = np.zeros((3,totpans))
    
    for i in range(0,totpans):
        forces[:,i] = Cp[i]*Areas[i]*np.divide(-normals[:,i],np.linalg.norm(normals[:,i]))
    
    
    return 
    
    
def MeshCentre(verts, surfs):
    
    return

    
    
verts, surfs = ml.cube(1)

normals = CalcNormals(verts, surfs)