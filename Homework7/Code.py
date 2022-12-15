# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 20:50:52 2022

@author: bene-
"""

import numpy as np

J=0.1
N=23
counter=0
cluster=np.zeros((N,N))
discoverd=np.full((N,N),False)
def int_lattice(N):
    return np.random.choice([-1, 1], size=(N, N))
    
def check_p(lattice):
    return np.where(np.random.rand()<(1-np.exp(2*J)*lattice),1,0)

def check_bonds(lattice):
    N=lattice.shape[0]
    
    rolled1=np.roll(lattice,1,axis=1)
    rolled2=np.transpose(np.roll(np.transpose(lattice),1,axis=1))
    
    bond1=check_p(np.abs(rolled1+lattice)/2)
    bond2=check_p(np.abs(rolled2+lattice)/2)
    
    return bond1,bond2

def ident_cluster(bond1,bond2):
    global counter
    global discoverd
    global cluster
    for i in range(N):
        for j in range(N):
            if (discoverd[i,j]==False):
                discoverd[i,j]=True
                cluster[i,j]=counter
                check_right(bond1,bond2,i,j)
                check_down(bond1,bond2,i,j)
                check_up(bond1,bond2,i,j)
                check_left(bond1,bond2,i,j)
                counter+=1
    
            
    
def check_right(bond1,bond2,i,j):
    global counter
    global cluster
    global discoverd
    if (bond1[i,j]==1): 
        cluster[(i+1)%N,j]=counter
        if (discoverd[(i+1)%N,j]==False):
            discoverd[(i+1)%N,j]=True
            check_right(bond1,bond2,(i+1)%N,j)
            check_down(bond1,bond2,(i+1)%N,j)
            check_up(bond1,bond2,(i+1)%N,j)
            
            
def check_left(bond1,bond2,i,j):
    global counter
    global cluster
    global discoverd
    if (bond1[(i-1)%N,j]==1): 
        cluster[(i-1)%N,j]=counter
        if (discoverd[(i-1)%N,j]==False):
            discoverd[(i-1)%N,j]=True
            check_down(bond1,bond2,(i-1)%N,j)
            check_up(bond1,bond2,(i-1)%N,j)
            check_left(bond1,bond2,(i-1)%N,j)

    
    
def check_down(bond1,bond2,i,j):
     global counter
     global cluster
     global discoverd
     if (bond2[i,j]==1):
         cluster[i,(j+1)%N]=counter
         if (discoverd[i,(j+1)%N]==False):
             discoverd[i,(j+1)%N]=True
             check_down(bond1,bond2,i,(j+1)%N)
             check_right(bond1,bond2,i,(j+1)%N)
             check_left(bond1,bond2,i,(j+1)%N)
             
    
def check_up(bond1,bond2,i,j):
     global counter
     global cluster
     global discoverd
     if (bond2[i,(j-1)%N]==1):
         cluster[i,(j-1)%N]=counter
         if (discoverd[i,(j-1)%N]==False):
             discoverd[i,(j-1)%N]=True
             check_down(bond1,bond2,i,(j-1)%N)
             check_right(bond1,bond2,i,(j-1)%N)
             check_left(bond1,bond2,i,(j-1)%N)

           
         
    

for i in range (20000):
    lattice=int_lattice(N)
    bond1,bond2=check_bonds(lattice)
   # print(bond1)
   # print(bond2)
    ident_cluster(bond1,bond2)
   # print(np.transpose(cluster))
    