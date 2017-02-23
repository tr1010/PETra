#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshes for Newtonian impact panel method
Created on Fri Jan 20 13:11:44 2017

@author: tr1010
import numpy as np
"""
import numpy as np

def CalcNormals(verts, surfs):
    """
    CalcNormals calculates the outward pointing unit normal vector for each of
    the quadrilateral surfaces making up the shape.
    
    Inputs:
        verts: a 3xn numpy array of reals describing the positions of the shape's 
               vertices 
        surfs: a 4xn numpy array of integers where each column gives the indices
               of the 4 vertices in verts variable that make up a surface of the
               shape. They are labelled in counterclockwise convention
    Outputs:
        normals: a 3xn numpy array of reals describing the unit normal vectors 
                 of each of the n surfaces making up the shape
    """
    
    totpans = np.size(surfs, 1)
    normals = np.zeros((3,totpans))
    
    for i in range(0,totpans):
        # Extract points making up the panel
        v1 = verts[:,surfs[1,i]] - verts[:,surfs[0,i]]
        v2 = verts[:,surfs[3,i]] - verts[:,surfs[0,i]]
        normals[:,i] = np.divide(np.cross(v1, v2),np.linalg.norm(np.cross(v1, v2)))
        
    return normals

def MeshCentres(verts, surfs):
    """
    MeshCentres calculates the centres of the quadrilateral surfaces which make
    make up the shape under consideration. The 'mesh centres' in this case are
    the points where the diagonal lines of the quadrilateral are nearest (in 3D
    space)
    
    Inputs:
        verts: a 3xn numpy array of reals describing the positions of the shape's 
               vertices 
        surfs: a 4xn numpy array of integers where each column gives the indices
               of the 4 vertices in verts variable that make up a surface of the
               shape. They are labelled in counterclockwise convention
               
    Outputs:
        centres: a 3xn array of reals describing the positions of each of the 
                 surfaces' centres
    
    """
    tot = np.size(surfs,1)
    centres = np.zeros((3,tot))
    bmps = np.zeros((3,4))
    for i in range(0,tot):
        #Find bimedian vectors
        bmps[:,0] = verts[:,surfs[0,i]] + 0.5*(verts[:,surfs[1,i]] - verts[:,surfs[0,i]])
        bmps[:,1] = verts[:,surfs[1,i]] + 0.5*(verts[:,surfs[2,i]] - verts[:,surfs[1,i]])
        bmps[:,2] = verts[:,surfs[2,i]] + 0.5*(verts[:,surfs[3,i]] - verts[:,surfs[2,i]])
        bmps[:,3] = verts[:,surfs[3,i]] + 0.5*(verts[:,surfs[0,i]] - verts[:,surfs[3,i]])
        
        A = np.transpose(np.array([bmps[:,2] - bmps[:,0],bmps[:,1] - bmps[:,3]]))
        b = np.transpose(np.subtract(bmps[:,1],bmps[:,0]))
        ts,resi,rank,s = np.linalg.lstsq(A,b)
        
        # Closest point on each line
        L1 = bmps[:,0] + (bmps[:,2] - bmps[:,0])*ts[0]
        L2 = bmps[:,1] + (bmps[:,3] - bmps[:,1])*ts[1]

        # Now find average
        centres[:,i] = np.divide(np.add(L1,L2),2)
        
    return centres

def MeshAreas(verts,surfs):
    """
    MeshAreas calculates the areas of the shape's surfaces. It assumes that 
    the surfaces are all quadrilaterals 
    
    Inputs:
        verts: a 3xn numpy array of reals describing the positions of the shape's 
               vertices 
        surfs: a 4xn numpy array of integers where each column gives the indices
               of the 4 vertices in verts variable that make up a surface of the
               shape. They are labelled in counterclockwise convention
    
    Outputs:
        areas: a n-element array of the areas of each of the shape's surfaces.
        
    """
    tot = np.size(surfs,1)
    areas = np.zeros((tot,1))
    
    for i in range(0,tot):
        # Find diagonal vectors
        p = verts[:,surfs[2,i]] - verts[:,surfs[0,i]]
        q = verts[:,surfs[3,i]] - verts[:,surfs[1,i]]
    
        areas[i] = 0.5*np.linalg.norm(np.cross(p,q)) 
    
    return areas

def Box(L,rho,CG_off=np.zeros(3)):
    """
    Box generates a Cuboid geometry for use in the the trajectory solver
    
    Inputs:  
        L : a 3-element array giving the x,y,z lengths of the cuboid      
        rho: density of the cuboid material (in kg/m3). Box assumes the cuboid 
             is a solid of uniform density       
        CG_off: Ignore this. Default a 3-element array of zeros
             
    Outputs:     
        verts: a 3x8 numpy array of reals describing the positions of the box's 
               vertices          
        surfs: a 4x6 numpy array of integers where each column gives the indices
               of the 4 vertices in verts variable that make up a surface of the
               cuboid. They are labelled in counterclockwise convention        
        areas: a 6-element array of the areas of each of the cuboid's surfaces
        normals: a 3x6 numpy array of reals. Each column describes the unit
                 normal vector (pointing outwards) for each of the cuboid's 
                 surfaces.            
        centroids: a 3x6 numpy array of reals describing the position of each 
                   surface's centroid.            
        CoG: Vector describing position of cuboid's centre of gravity. By 
             default is always zeros (coordinate system origin is at CoG)]       
        I: 3x3 numpy array of reals describing the cuboid's inertia tensor
        mass:real describing the cuboid's mass
    """
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