# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 13:15:38 2022

@author: bene-
"""

import numpy as np
from  math import *
import matplotlib.pyplot as plt
import scipy.optimize as so


N=20
J=1
T=1
h=0.2

def init_1D(N):
    lattice=np.zeros(N)
    return lattice

def step(lattice,N): #Does one stop in the Markov-Chain and return the next constallation 
    global J
    global T           
    beta=1/T 
    M=np.random.randint(0,N) #selects a random sin in the cain 
    if np.random.rand()>0.5:#generates a new random orientation of the spin 
        spin_new=1
    else:
        spin_new=-1 
    delta_H=calc_delta_H(lattice,M,spin_new,N,h) #calculates the energy difference from the old an new constallation 
    if delta_H<0:
        lattice[M]=spin_new
    else:
        if np.exp(-beta * delta_H)>np.random.rand():
            lattice[M]=spin_new
    return(lattice)
    

def calc_delta_H(lattice,M,spin_new,N,h):#calculates the energy difference from the old an new constallation 
    global J 
    spin_old=lattice[M]
    delta_H=0
    delta_H+= -J*(lattice[(M-1)%N]*spin_old - lattice[(M-1)%N]*spin_new)
    delta_H+= -J*(lattice[(M+1)%N]*spin_old - lattice[(M+1)%N]*spin_new)   
    delta_H+=-h*(spin_old-spin_new) 
    return -delta_H

def calc_magn (lattice,N): #calculates the magentization for one canstallation
    M=0
    M=np.sum(lattice)
    return M/N

def show_lattice(lattice,N): #visialization of the lattice
    global T
    global J 
    
    X ,Y= np.meshgrid(np.arange(N),1)
    U = np.cos(0.5*np.pi*lattice)
    V = np.sin(0.5*np.pi*lattice)
    plt.figure(figsize=(6,6), dpi=400)
    Q = plt.quiver(X, Y, U, V, units='width')
    plt.quiverkey(Q, 0.1, 0.1, 1, r'.' , labelpos='E',
                    coordinates='figure')
    
    plt.axis('off')
 #   plt.savefig('Abbildungen/Gitter_%s.png' %(str(T)))
    plt.show()
    
def theo_mag(N,h):
    global T
    
    return (T/N)*(N*np.exp((J/T))*(np.sinh(h/T)/T - (np.sinh(h/T)*np.cosh(h/T))/(T*np.sqrt(np.exp((-(4*J)/T)) + np.sinh(h/T)**2)))*(np.exp((J/T))*(np.cosh(h/T) - np.sqrt(np.exp((-(4*J)/T)) + np.sinh(h/T)**2)))**(N - 1)+N*np.exp((J/T))*(np.sinh(h/T)/T + (np.sinh(h/T)*np.cosh(h/T))/(T*np.sqrt(np.exp((-(4*J)/T)) + np.sinh(h/T)**2)))*(np.exp((J/T))*(np.cosh(h/T) + np.sqrt(np.exp((-(4*J)/T)) + np.sinh(h/T)**2)))**(N - 1))/((np.exp((J/T))*(np.cosh(h/T)+np.sqrt(np.sinh(h/T)**2+np.exp((-4*(J/T))))))**N+(np.exp((J/T))*(np.cosh(h/T)-np.sqrt(np.sinh(h/T)**2+np.exp((-4*(J/T))))))**N)


E=20 #amount of points for the graph
I=200 #initialization
A=100  #autocorrelationcorrection
C=5000 #Number of components for an ensemble
D=5 #Number of runs for error estimate

final1error=np.zeros(D*E)
final2error=np.zeros(D*(E+1))
N1=np.arange(1,21,1)
theo_magN=np.zeros(E)#calcultaes the theo_mag for every pont in the graph s.t. it doesnt have to be calculated every time when used
for i in range(E):
    theo_magN[i]=theo_mag(N1[i],h)
plt.figure()
plt.plot(N1,theo_magN,label='<m>_{theo}',marker='o')
plt.xlabel(r"$N$")
plt.ylabel(r"$m(N)$")


final_M=np.zeros(E)
M=0
colours=['b','g','r','c','m','y','k']

for m in range(D):
    print(m)
    for j in range(E): 
        N=j+1
        lattice=init_1D(N)
        for i in range(I): 
            (lattice)=step(lattice,N)
        for k in range(C):
            for i in range(A):
                (lattice)=step(lattice,N)
            M+=calc_magn(lattice,N)
        print(N)
        final_M[j]=M/C
        #show_lattice(lattice,N)
        M=0
        popt, pcov=so.curve_fit(theo_mag,N1,final_M,maxfev=2000)
        perr = np.sqrt(np.diag(pcov))
        print(perr)
        final1error[j+m*E]=(final_M[j]-theo_magN[j])**2 #saving the error for every calculated magnetization
    plt.plot(N1,final_M,color=colours[m+1],marker='o')
sigmaN=np.sqrt(1/(D*E)*np.sum(final1error)) #calculates the sigma of from the errors saved above 
print(sigmaN)
plt.legend()
plt.show() 
  
h1=np.arange(-1.0,1.1,0.1)
theo_magh=np.zeros(E+1)
for i in range(E+1):
    theo_magh[i]=theo_mag(N,h1[i])
plt.figure()
plt.plot(h1,theo_magh,label='<m>_{theo}',marker='o')
plt.xlabel(r"$h$")
plt.ylabel(r"$m(h)$")
    
N=20    

final_M=np.zeros(E+1)
for m in range(D):
    print(m)
    for j in range(E+1):
        h=-1+(j/10)
        lattice=init_1D(N)
        for i in range(I): 
            (lattice)=step(lattice,N)
        for k in range(C):
            for i in range(A):
                (lattice)=step(lattice,N)
            M+=calc_magn(lattice,N)
        print(h)
        final_M[j]=M/C
        #show_lattice(lattice,N)
        M=0
        final2error[j+m*(E+1)]=(final_M[j]-theo_magh[j])**2
    plt.plot(h1,final_M,color=colours[m+1],marker='o')
plt.legend()        
plt.show()  
sigmah=np.sqrt(1/(D*(E+1))*np.sum(final2error))
print(sigmah)     
