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
"""
import numpy as np

# cube
def cube(L):
    verts = L*np.array([[0., 1., 1., 0., 0., 1., 1., 0.],
                        [0., 0., 0., 0., 1., 1., 1., 1.],
                        [0., 0., 1., 1., 0., 0., 1., 1.]])
    
    surfs = np.array([[0, 1, 1, 0, 2, 6],
                      [1, 5, 0, 3, 6, 5],
                      [2, 6, 4, 7, 7, 4],
                      [3, 2, 5, 4, 3, 7]])
    
    return verts, surfs
