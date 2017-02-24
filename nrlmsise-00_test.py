#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nrlmsise-00 test
Created on Fri Feb 24 11:42:55 2017

@author: tr1010
"""

import sys  
sys.path.append('NRLMSISE-00/Python-NRLMSISE-00-master')
import numpy as np
from nrlmsise_00_header import *
from nrlmsise_00 import *

def gtd7_test():
    output = [nrlmsise_output() for _ in range(17)]
    Input = [nrlmsise_input() for _ in range(17)]
    flags = nrlmsise_flags()
    aph = ap_array()
    
    for i in range(7):
        aph.a[i]=100
    
    flags.switches[0] = 1 #output in m rather than cm

    for i in range(1,24):
        flags.switches[i]=1

    for i in range(17):
        Input[i].doy=172
        Input[i].year=0; #/* without effect */
        Input[i].sec=29000;
        Input[i].alt=400;
        Input[i].g_lat=60;
        Input[i].g_long=-70;
        Input[i].lst=16;
        Input[i].f107A=150;
        Input[i].f107=150;
        Input[i].ap=4;
	
    Input[1].doy=81;
    Input[2].sec=75000;
    Input[2].alt=1000;
    Input[3].alt=100;
    Input[10].alt=0;
    Input[11].alt=10;
    Input[12].alt=30;
    Input[13].alt=50;
    Input[14].alt=70;
    Input[16].alt=100;
    Input[4].g_lat=0;
    Input[5].g_long=0;
    Input[6].lst=4;
    Input[7].f107A=70;
    Input[8].f107=180;
    Input[9].ap=40;

    Input[15].ap_a = aph
    Input[16].ap_a = aph

    #evaluate 0 to 14
    for i in range(15):
        gtd7(Input[i], flags, output[i])

    #/* evaluate 15 and 16 */
    flags.switches[9] = -1
    for i in range(15, 17):
        gtd7(Input[i], flags, output[i])
    
    #/* output type 2 */
    for i in range(3):
        print('\n', end='')
        print("\nDAY   ", end='')
        for j in range(5):
            print("         %3i" % Input[i*5+j].doy, end='')
        print("\nUT    ", end='')
        for j in range(5):
            print("       %5.0f" % Input[i*5+j].sec, end='')
        print("\nALT   ", end='')
        for j in range(5):
            print("        %4.0f" % Input[i*5+j].alt, end='')
        print("\nLAT   ", end='')
        for j in range(5):
            print("         %3.0f" % Input[i*5+j].g_lat, end='')
        print("\nLONG  ", end='')
        for j in range(5):
            print("         %3.0f" % Input[i*5+j].g_long, end='')
        print("\nLST   ", end='')
        for j in range(5):
            print("       %5.0f" % Input[i*5+j].lst, end='')
        print("\nF107A ", end='')
        for j in range(5):
            print("         %3.0f" % Input[i*5+j].f107A, end='')
        print("\nF107  ", end='')
        for j in range(5):
            print("         %3.0f" % Input[i*5+j].f107, end='')

        print('\n\n', end='')
        
        print("\nTINF  ", end='')
        for j in range(5):
            print("     %7.2f" % output[i*5+j].t[0], end='')
        print("\nTG    ", end='')
        for j in range(5):
            print("     %7.2f" % output[i*5+j].t[1], end='')
        print("\nHE    ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[0], end='')
        print("\nO     ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[1], end='')
        print("\nN2    ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[2], end='')
        print("\nO2    ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[3], end='')
        print("\nAR    ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[4], end='')
        print("\nH     ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[6], end='')
        print("\nN     ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[7], end='')
        print("\nANM   ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[8], end='')
        print("\nRHO   ", end='')
        for j in range(5):
            print("   %1.3e" % output[i*5+j].d[5], end='')
        print('\n')
    
    return
    
gtd7_test()