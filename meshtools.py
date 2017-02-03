#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshes for Newtonian impact panel method
Created on Fri Jan 20 13:11:44 2017

@author: tr1010

Here are some meshes for a number of canonical shapes which can then be used
for Newtonian impact theory analyses
Pyramid thing
Sphere
Cylinder
Box
Flat plate
These correspond to the objects used by ESA's DRAMA-SESAM Object-oriented 
demise tool


TO DO:
    - Find a way to distinguish between 'internal' and 'external' panels (in order to model things like e.g. solar panels)
    - Implement a coordinate system origin for each mesh (to track orientation of object)
"""
import numpy as np

# Function for calculating normal vectors (anticlockwise numbering convention)
def CalcNormals(verts, surfs):
    
    totpans = np.size(surfs, 1)
    normals = np.zeros((3,totpans))
    
    for i in range(0,totpans):
        # Extract points making up the panel
        v1 = verts[:,surfs[1,i]] - verts[:,surfs[0,i]]
        v2 = verts[:,surfs[3,i]] - verts[:,surfs[0,i]]
        normals[:,i] = np.divide(np.cross(v1, v2),np.linalg.norm(np.cross(v1, v2)))
        
    return normals

# Function for calculating geometric centres of meshes using trapezium relations
def MeshCentres(Verts, surfs):
    tot = np.size(surfs,1)
    Centres = np.zeros((3,tot))
    bmps = np.zeros((3,4))
    for i in range(0,tot):
        #Find bimedian vectors
        bmps[:,0] = Verts[:,surfs[0,i]] + 0.5*(Verts[:,surfs[1,i]] - Verts[:,surfs[0,i]])
        bmps[:,1] = Verts[:,surfs[1,i]] + 0.5*(Verts[:,surfs[2,i]] - Verts[:,surfs[1,i]])
        bmps[:,2] = Verts[:,surfs[2,i]] + 0.5*(Verts[:,surfs[3,i]] - Verts[:,surfs[2,i]])
        bmps[:,3] = Verts[:,surfs[3,i]] + 0.5*(Verts[:,surfs[0,i]] - Verts[:,surfs[3,i]])
        
        A = np.transpose(np.array([bmps[:,2] - bmps[:,0],bmps[:,1] - bmps[:,3]]))
        b = np.transpose(np.subtract(bmps[:,1],bmps[:,0]))
        ts,resi,rank,s = np.linalg.lstsq(A,b)
        
        # Closest point on each line
        L1 = bmps[:,0] + (bmps[:,2] - bmps[:,0])*ts[0]
        L2 = bmps[:,1] + (bmps[:,3] - bmps[:,1])*ts[1]

        # Now find average
        Centres[:,i] = np.divide(np.add(L1,L2),2)
        
    return Centres

# Function for calculating Areas of meshes using trapezium relations
def MeshAreas(Verts,surfs):
    tot = np.size(surfs,1)
    Areas = np.zeros((tot,1))
    
    for i in range(0,tot):
        # Find diagonal vectors
        p = Verts[:,surfs[2,i]] - Verts[:,surfs[0,i]]
        q = Verts[:,surfs[3,i]] - Verts[:,surfs[1,i]]
    
        Areas[i] = 0.5*np.linalg.norm(np.cross(p,q)) 
    
    return Areas

# cube
def Box(L,rho,CG_off):
    verts = np.array([[0., 1.*L[0], 1.*L[0], 0., 0., 1.*L[0], 1.*L[0], 0.],
                      [0., 0., 0., 0., 1.*L[1], 1.*L[1], 1.*L[1], 1.*L[1]],
                      [0., 0., 1.*L[2], 1.*L[2], 0., 0., 1.*L[2], 1.*L[2]]])
    
    surfs = np.array([[0, 1, 1, 0, 2, 6],
                      [1, 5, 0, 3, 6, 5],
                      [2, 6, 4, 7, 7, 4],
                      [3, 2, 5, 4, 3, 7]])
    
    #Centre of gravity with offset given by user
    CoG = np.add(L/2,CG_off)
    
    # Move origin to centre of gravity
    for i in range(0,np.size(verts,1)):
        verts[:,i] = np.subtract(verts[:,i],CoG)
    
    # New CoG is at origin
    CoG = np.zeros((3,1))
    
    # Now calculate inertia tensor
    # for a cuboid:
    mass = rho*np.prod(L)
    temp = mass/12.
    I = np.array([[temp*(L[1]**2 + L[2]**2), 0., 0.],
                  [0., temp*(L[0]**2 + L[2]**2), 0.],
                  [0., 0., temp*(L[0]**2 + L[1]**2)]])
    
    # Calculate new principle axis at offset centre of mass
    I = I + mass*(np.eye(3)*np.linalg.norm(CG_off)**2 - np.outer(CG_off,CG_off))
    
    normals = CalcNormals(verts,surfs)
    
    areas = MeshAreas(verts,surfs)
    
    centroids = MeshCentres(verts,surfs)
    
    
    return verts, surfs, areas, normals, centroids, CoG, I, mass