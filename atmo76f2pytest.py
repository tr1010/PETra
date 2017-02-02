#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:16:58 2017

@author: tr1010
"""

import atmospy76 as atm
import numpy as np
sigma = np.float64()
delta = np.float64()
theta = np.float64()
alt = 200.
sigma,delta,theta = atm.atmosphere(alt)

