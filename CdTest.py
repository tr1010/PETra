#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cd Test: test script for Cd_Calc function
Created on Thu Dec  1 11:37:27 2016

@author: tr1010
"""
import Trajectory as tr

M = 8.
Kn = 0.001
V = 3000
RT = 40000

Cd = tr.Cd_Calc(M,Kn,V,RT)