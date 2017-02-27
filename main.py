#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PETra -- main function

Created on Tue Nov 29 18:05:23 2016

@author: tr1010 (Thomas Rees)
"""
import scipy as sp
import trajectory as tr
import postprocessing as pp
import inputs as InFile
import input_processing as PI
    
def main():
    """
    PETra's main function calls three functions:
        process_inputs
        sp.integrate.odeint (a RK45 ODE integrator included in scipy)]
        postprocess
    It returns two Dictionaries: InputDict and ResultsDict which contain all
    inputs used in the run and all the outputs generated
    """
    
    earth, mass, areas, normals, centroids, I, t, x0, scLengthScale, aero_params, InputDict = PI.process_inputs(InFile.get_inputs)

    sol = sp.integrate.odeint(tr.traj_uvw,x0,t,args=(earth, mass, areas, normals, centroids, I, scLengthScale, aero_params))
    
    ResultsDict = pp.postprocess(t, sol, earth, mass, areas, normals, centroids, I, x0, scLengthScale)
    
    return InputDict, ResultsDict
    
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    InputDict, ResultsDict = main()   
    
