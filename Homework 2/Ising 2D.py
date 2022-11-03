# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 12:36:03 2022

@author: bene-
"""

import numpy as np
from  math import *
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import mpmath as mp

N_x=10
N_y=10
h=0.0
T=1
J=0.5

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
                if delta_H<0:
                    lattice[M,N]=spin_new
                else:
                    if np.exp(-beta * delta_H)>np.random.rand():
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
    
def show_lattice(lattice): #gibt die Matrix des Gitters als Pfeile auf den Gitterpunkten aus. 
    global N_x
    global N_y
    global J

    X, Y = np.meshgrid(np.arange(0,N_x),np.arange(0,N_y))
    U = np.cos(0.5*np.pi*lattice)
    V = np.sin(0.5*np.pi*lattice)
    plt.figure(figsize=(6,6), dpi=400)
    Q = plt.quiver(X, Y, U, V, units='width')
    plt.quiverkey(Q, 0.1, 0.1, 1, r'.' , labelpos='E',
                    coordinates='figure')
    
    plt.axis('off')
#    plt.savefig('Abbildungen/Gitter_%s.png' %(str(T)))
    plt.show()

def calc_magn (lattice): #Berechnet die Magnetisierung eines Gitters
    global N_x
    global N_y

    M=0
    M=np.sum(lattice)
    return M/(N_x*N_y)

def calc_eps(lattice):
    global N_x
    global N_y
    global h
    global J
    
    E=0.0
    rolled1=np.roll(lattice,1,axis=1)
    rolled2=np.roll(lattice,-1,axis=1)#verschiebt Zeilen R weit nach rechts bzw. links
    rolled3=np.transpose(np.roll(np.transpose(lattice),1,axis=1))
    rolled4=np.transpose(np.roll(np.transpose(lattice),-1,axis=1))    
    E+=-J*np.sum(np.multiply(rolled1,lattice))
    E+=-J*np.sum(np.multiply(rolled2,lattice))
    E+=-J*np.sum(np.multiply(rolled3,lattice))
    E+=-J*np.sum(np.multiply(rolled4,lattice))
    E+=-h*np.sum(lattice)

    return E/(N_x*N_y)

def eps_theo():
    global J

    A=mp.sech(2*J)**2
    B=mp.tanh(2*J)**2

    return -2*J*mp.coth(2*J)*(1+(2/np.pi)*(2*(mp.tanh(2*J)**2)-1)*mp.ellipk(4*A*B))
    
def abs_m_theo():
    global J
    if J>0.4406867935:
        return (1-(1/mp.sinh(2*J)**4))**(1/8)
    else: 
        return 0
    
def m_h_5():
    global N_x
    global N_y
    global h
    E=20 #amount of points for the graph
    I=10 #initialization
    A=10 #autocorrelationcorrection
    C=500 #Number of components for an ensemble
    D=5
    final_M=np.zeros(E+1) 
    M=0
    lattice=init_2D()
    colours=['b','g','r','c','m','y','k']    
    for m in range(D):       
        N_x=N_y=(m+1)*4
        final_M=np.zeros(E+1)
        for j in range(E+1):
            h=-1+(j/10)
            lattice=init_2D()
            for i in range(I): 
                (lattice)=step(lattice)
            for k in range(C):
                for i in range(A):
                    (lattice)=step(lattice)
                M+=calc_magn(lattice)
            print(h)
            final_M[j]=M/C
            M=0
            #show_lattice(lattice)
        plt.plot(np.arange(-1.0,1.1,0.1),final_M,color=colours[m+1],label='$N_x=N_y=$%s'%str((m+1)*4))
    plt.legend()
    plt.xlabel(r"$h$")
    plt.ylabel(r"$m(h)$")
    plt.show()
    
def eps_J_6():
    global N_x
    global N_y
    global h
    global J
    

    E=16 #amount of points for the graph
    I=10 #initialization
    A=10 #autocorrelationcorrection
    C=100 #Number of components for an ensemble
    D=5
    lattice=init_2D()
    colours=['b','g','r','c','m','y','k']   
    h=0.0
    final_eps=np.zeros(E)
    eps_theo_array=np.zeros(E)
    for i in range(E):
        J=0.25+(i*0.75/(E-1))
        eps_theo_array[i]=eps_theo()
        
        
    for m in range(D):       
        N_x=N_y=(m+1)*4
        final_eps=np.zeros(E)
        for j in range(E):
            J=0.25+(j*0.75/(E-1))
            lattice=init_2D()
            En=0
            for i in range(I): 
                (lattice)=step(lattice)
            for k in range(C):
                for i in range(A):
                    (lattice)=step(lattice)
                En+=calc_eps(lattice)
            print(J)
            final_eps[j]=En/C
            #show_lattice(lattice)
        plt.plot(np.arange(0.25,1.05,0.05),final_eps,color=colours[m+1],label='$N_x=N_y=$%s'%str((m+1)*4))
    plt.plot(np.arange(0.25,1.05,0.05),eps_theo_array,label='$N_x=N_y=\infty$ (theo)')
    plt.legend()
    plt.xlabel(r"$J$")
    plt.ylabel(r"$\epsilon (J)$")
    plt.show()  
    
    
def abs_m_J_7():
    global N_x
    global N_y
    global h
    global J
    

    E=16 #amount of points for the graph
    I=10 #initialization
    A=10 #autocorrelationcorrection
    C=1000 #Number of components for an ensemble
    D=5
    lattice=init_2D()
    colours=['b','g','r','c','m','y','k']   
    h=0.0
    final=np.zeros(E)
    abs_m_theo_array=np.zeros(E)
    for i in range(E):
        J=0.25+(i*0.75/(E-1))
        abs_m_theo_array[i]=abs_m_theo()
        
        
    for m in range(D):       
        N_x=N_y=(m+1)*4
        final=np.zeros(E)
        for j in range(E):
            J=0.25+(j*0.75/(E-1))
            lattice=init_2D()
            abs_m=0
            for i in range(I): 
                (lattice)=step(lattice)
            for k in range(C):
                for i in range(A):
                    (lattice)=step(lattice)
                abs_m+=np.abs(calc_magn(lattice))
            print(J)
            final[j]=abs_m/C
            #show_lattice(lattice)
        plt.plot(np.arange(0.25,1.05,0.05)**(-1),final,color=colours[m+1],label='$N_x=N_y=$%s'%str((m+1)*4))
    print(abs_m_theo_array)
    plt.plot(np.arange(0.25,1.05,0.05)**(-1),abs_m_theo_array,label='$N_x=N_y=\infty$ (theo)')
    # plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.legend()
    plt.xlabel(r"$J^{-1}$")
    plt.ylabel(r"$\langle |m|\rangle (J)$")
    plt.show()  
   
# m_h_5()
# eps_J_6()
abs_m_J_7()
        
    
