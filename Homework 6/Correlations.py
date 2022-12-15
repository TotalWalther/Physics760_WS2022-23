# -*- coding: utf-8 -*-



import numpy as np
from  math import *
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import mpmath as mp
import random


N_x=23
N_y=23
h=0.0
T=1
J=0.43/T



def init_2D(): #erzeugt ein LxL-Gitter, wobei ein initialer Kaltstart gewählt wird 
    global N_x
    global N_y
    lattice = np.zeros((N_x,N_y))
    return lattice


def step(lattice): #Does one stop in the Markov-Chain and return the next constallation 
    global N_x
    global N_y
    global J
    global T           
    beta=1/T 
    for i in range (N_x):
        M=i
        for j in range (N_y):
                N=j #selects a random sin in the cain 
                if np.random.rand()>0.5:#generates a new random orientation of the spin 
                    spin_new=1
                else:
                    spin_new=-1 
                delta_H=calc_delta_H(lattice,M,N,spin_new) #calculates the energy difference from the old an new constallation 
                P_acc = np.exp(float(-delta_H))
                if P_acc>np.random.rand():
                    lattice[M,N]=spin_new
    return(lattice)

def calc_delta_H(lattice,M,N,spin_new):
    global J
    global N_x
    global N_y
    global h

    spin_old=lattice[M][N]
    delta_H=0
    delta_H+=J* (lattice[(M-1)%(N_x)][N] * spin_new) \
                    -J* (lattice[(M-1)%(N_x)][N] * spin_old)
    delta_H+=J* (lattice[(M+1)%(N_x)][N] * spin_new) \
                    -J* (lattice[(M+1)%(N_x)][N] * spin_old)  
    delta_H+=J* (lattice[M][(N-1)%(N_y)] * spin_new) \
                    -J* (lattice[M][(N-1)%(N_y)] * spin_old)   
    delta_H+=J* (lattice[M][(N+1)%(N_y)] * spin_new) \
                    -J* (lattice[M][(N+1)%(N_y)] * spin_old) 
    delta_H+=h*(spin_new-spin_old)
    #Berechnet jeweils den Energieunterschied für benachbarte Spins für den alten und neuen Spin
    return -delta_H

def calc_C(lattice):
    lattice_k=np.fft.fft2(lattice, norm="ortho")
    lattice_mk=np.fft.ifft2(lattice, norm="ortho")
    C=1/N_x*np.fft.ifft2(lattice_k*lattice_mk, norm="ortho")
    return C

def calc_C_r(C,r):
    C_r=0
    for i in range(int(np.sqrt(r)+0.5)):
        j=int(np.sqrt(r-i**2)+.5)
        C_r+=C[i,j]+C[j,i]
        
    return C_r/(int(np.sqrt(r)+0.5)*2)
        
    
            
    
    
final_C_r=[]  
lattice=init_2D()
for j in range (1):
    for i in range(2000):
        lattice=step(lattice)
    C=np.real(calc_C(lattice))
    print(np.real(calc_C_r(C,5)))
for i in range(N_x):
    if i==0: continue
    final_C_r.append(np.real(calc_C_r(C,i)))
plt.plot(final_C_r)
#    plt.plot(np.real(C_r[0]))
#     final_C_r.append(calc_C_r(lattice,0))   
# print(np.average(final_C_r))
    