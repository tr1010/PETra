#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trajectories main
Created on Tue Nov 29 18:05:23 2016

@author: tr1010
"""
import scipy as sp
import Trajectory as tr
import postprocessing as pp
import InputFile as InFile
import Process_Inputs as PI
    
def main(): 
    
    earth, mass, areas, normals, centroids, I, t, x0, scLengthScale, InputDict = PI.process_inputs(InFile.get_inputs)

    sol = sp.integrate.odeint(tr.traj_uvw,x0,t,args=(earth, mass, areas, normals, centroids, I, scLengthScale))
    
    ResultsDict = pp.postprocess(t, sol, earth, mass, areas, normals, centroids, I, x0, scLengthScale)
    
    return InputDict, ResultsDict
    
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    InputDict, ResultsDict = main()   
    
