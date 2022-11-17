# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 14:30:05 2022

@author: bene-
"""

import numpy as np
from  math import *
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import mpmath as mp


beta=1
N=20
h=0.5
J=1/N

def leapfrog(p_0,phi_0,N_md):
    global beta 
    global J
    global N 
    global h
    
    eps=1/N_md
    phi=phi_0+eps/2*p_0
    p=p_0
    for i in range(N_md-1):
        p=p-eps*(phi/(beta*J)-N*mp.tanh(beta*h+phi))
        phi=phi+eps*p
    p=p-eps*(phi/(beta*J)-N*mp.tanh(beta*h+phi))
    phi=phi+eps/2*p
    
    return [p,phi]

def H(p,phi):
    global beta 
    global J 
    global h 
    
    return p**2/2+phi**2/(2*beta*J)-N*mp.log(2*mp.cosh(beta*h+phi))
    
p_0=0.5
phi_0=.5
H_0=H(p_0,phi_0)
for i in range(100):
    p,phi=leapfrog(p_0,phi_0,i+2)
    plt.plot(i+2,abs((H(p,phi)-H_0)/H_0), 'x', color='b')
plt.semilogy()    
plt.show()    