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
h=-1

def init_1D():
    global N
    lattice=np.zeros(N)
    return lattice

def step(lattice): #führt einen Schritt in der Markov-Kette aus und gibt das nächste Gitter der Kette zurück
    global N
    global T
    global J
    beta=1/T 
    M=np.random.randint(0,N) #wählt einen zufälligen Gitterpunkt aus
    if np.random.rand()>0.5:
        spin_new=1
    else:
        spin_new=-1 #zufällige neue Ausrichtung des Spins
    delta_H=calc_delta_H(lattice,M,spin_new) #Berechnet den Energieunterschied der alten und neuen Konfiguration
    if delta_H<0:
        lattice[M]=spin_new
    else:
        if np.exp(-beta * delta_H)>np.random.rand():
            lattice[M]=spin_new
    return(lattice)
    

def calc_delta_H(lattice,M,spin_new):
    global N
    global J 
    global h
    spin_old=lattice[M]
    delta_H=0
    delta_H+= -J*(lattice[(M-1)%N]*spin_new - lattice[(M-1)%N]*spin_old)
    delta_H+= -J*(lattice[(M+1)%N]*spin_old - lattice[(M+1)%N]*spin_old)   
    delta_H+=-h*(spin_new-spin_old) 

    #Berechnet jeweils den Energieunterschied für benachbarte Spins für den alten und neuen Spin
    return delta_H

def calc_magn (lattice): #Berechnet die Magnetisierung eines Gitters
    global N
    M=0
    M=np.sum(lattice)
    return M/N

def show_lattice(lattice): #gibt die Matrix des Gitters als Pfeile auf den Gitterpunkten aus. 
    global N 
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

E=40 #anzahl der temperaturen
I=50000 #initilisierung
A=100  #autokorrelationkorekktur
C=10    #temperaturintervall
D=30    #temperaturintervall

lattice=init_1D()
for i in range(1000):
    lattice=step(lattice)
show_lattice(lattice)
final_M=np.zeros(E)
M=0
colours=['b','g','r','c','m','y','k']

# for j in range(E): #Schleife für die Anzahl der zu untersuchenden Temperaturen
#     global T
#     beta=1/((j+C)/D) #T=(j+C)/D
#     for i in range(I): #Anpassung des Gitters an die neue Temperatur
#         (lattice)=step(lattice)
#     for k in range(500):#Schleife für die Anzahl der verwendeten Ensemble pro Temperatur
#         for i in range(A):
#             (lattice)=step(lattice)
#         M+=calc_magn(lattice)
#     print(1/beta)
#     final_M[j]=M/500
#     #show_lattice(gitter)
#     M=0
T=(np.arange(E)+C)/D
plt.plot(T,final_M)


