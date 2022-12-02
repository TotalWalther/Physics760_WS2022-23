# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 08:27:58 2022

@author: bene-
"""

#!/usr/bin/env python
# coding: utf-8

import numpy as np
from  math import *
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import mpmath as mp
import random




# The graph looks like we expected it. It shows that for large $N_{md}$ the energy of the old and new configuration is almost the same. This will later be important to achieve a high acceptance probability.  

# In[35]:

N_md = 10#Leapfrog integration steps
N_cfg = 100000
beta=1000
phi=0
p_0 = 0 
ar=0



# This is just defining of variables. 

# In[36]:


def leapfrog_plot():
    global N_md
    p=0.2
    phi=0
    for i in range(100):
        H_0=H(p,phi)
        N_md=i*10+10
        p_new,phi_new=leapfrog(p,phi)
        plt.plot(i*10+10,abs((H(p_new,phi_new)-H_0)/H_0), 'x', color='b')
        #print(p,phi)
    plt.semilogy()    
    plt.show() 




# In[42]:


def leapfrog(p_l,phi_l):
    # p_0_1,p_0_2,p_0_3,
    global beta 
    global f
    global x
    global J
    global N 
    global h
    global N_md
    
    eps=.1/N_md
    phi_l=phi_l+eps/2*p_l
 
    for i in range(N_md-1):
        p_l=p_l-eps*(2*phi_l)*(1+1/(2+phi_l**2))
        phi_l=phi_l+eps*p_l
        #print(p)
    p_l=p_l-eps*(2*phi_l)*(1+1/(2+phi_l**2))
    phi_l=phi_l+eps/2*p_l
    #print(p)
    return p_l,phi_l


# input: p_0, phi_0; output: p_f,phi_f  
# code as explained on the sheet

# In[43]:


def H(p_h,phi_h):
    global beta 
    global J 
    global h 
    global f
    global x
    global delta_f
    #print(phi_h[0]+x*phi_h[1]+x**2*phi_h[2])
    return  p_h**2/2+phi_h**2+np.log(2+phi_h**2)


# input p,phi: ; output: H(p,phi)

# In[44]:


def HMC(): #Does one iteration of the Markov-Chain and return phi
    global N_md
    global p_0
    global phi
    #global p
    global ar
    
    p_h=np.random.normal(loc=00.0, scale=1.0)
    p_l,phi_l = leapfrog(p_h,phi)    
    
    
    P_acc = np.exp(float(H(p_h,phi)-H(p_l,phi_l)))
    
        
    if P_acc > np.random.rand(): 
        
        phi = phi_l
        ar=ar+1

   
def autocorrelation(o_chain,t,mu_o):
    c_t=0
    for i in range(N_cfg-t):
        c_t+=(o_chain[i]-mu_o)*(o_chain[i+t]-mu_o)
    c_t=c_t/(N_cfg-t)
    return c_t


def normalized_autocorrelation(o_chain,t):
    mu_o=np.average(o_chain)
    Gamma_t=autocorrelation(o_chain,t,mu_o)/autocorrelation(o_chain,0,mu_o)
    return Gamma_t

def int_auto_corr(o_chain):
    tau_int=0
    tau_int=0.5*normalized_autocorrelation(o_chain,0)
    for i in range(N_cfg):
        temp=normalized_autocorrelation(o_chain,i+1)
        if temp>0:
            tau_int+=temp
        else: return tau_int
    return tau_int
# Classical HMC-Algo, which returns the next element of the markov chain. Candidates are created with the leapfrog-algo. In our case it also keeps track of the acceptance probability with ar.

# In[45]:




def markov_chain():
    global N_cfg
    global x
    global phi
    global ar
    
    
    phi_chain=[]
    for i in range(N_cfg):
        HMC()
        phi_chain.append(phi)
        #print(ar/(i+1))
    print(ar/(N_cfg),'Acceptence rate')
    ar=0
    phi_chain=np.array(phi_chain)
    # plt.plot(phi_chain)
    o_chain=np.cos(np.sqrt(1+phi_chain**2))
    # plt.plot(o_chain)
    # print(np.average(o_chain))
    return o_chain

def binning(o_chain,tau_int):
    o_chain_new=[]
    tau_int=int(tau_int)
    o_chain=o_chain[0:N_cfg-N_cfg%tau_int]
    for i in range(int((N_cfg-N_cfg%tau_int)/tau_int)):
        o_chain_new.append(np.average(o_chain[i*tau_int:(i+1)*tau_int]))
    return o_chain_new
        
def print_nor_auto(o_chain,tau_int,label1,label2):
    N_cfg_temp=len(o_chain)
    normal_auto= []
    exp= []
    for i in range(N_cfg_temp):
        #print(i)
        #print(autocorrelation(o_chain,i,np.average(o_chain)))
        normal_auto.append(normalized_autocorrelation(o_chain,i))
        exp.append(np.exp(-i/tau_int))
    plt.plot(normal_auto,label=label1)
    plt.plot(exp,label=label2)
    plt.legend()
    plt.show()
    
def show_auto_dec(o_chain,tau_int):
    global N_cfg
    tau_int_final=[]
    tau_int_final.append(tau_int)
    for i in range(2,int(tau_int)):
        N_cfg=len(o_chain)
        o_chain_new=binning(o_chain,i)
        N_cfg=len(o_chain_new)
        #print(N_cfg,'\n %s'%str(i))
        tau_int_new=int_auto_corr(o_chain_new)
        #print(tau_int_new,'tau_int_after','\n %s'%str(i))
        #print_nor_auto(o_chain_new,tau_int,str(1),str(2))
        tau_int_final.append(tau_int_new)
    N_cfg=len(o_chain)    
    plt.plot(tau_int_final)
    plt.xlabel(r"$N_{bin}$")
    plt.ylabel(r"$\tau_{int}(N_{bin}$")
    plt.show()
    
def show_error_inc(o_chain,tau_int):
    boot_var_final=[]
    for i in range(2,3*int(tau_int)):
        o_chain_new=binning(o_chain,i)
        o_chain_boot=booti(o_chain_new, 200)
        boot_var_final.append(np.sqrt(np.var(o_chain_boot)))   
    plt.plot(boot_var_final)
    plt.xlabel(r"$N_{bin}$")
    plt.ylabel(r"$\sqrt{var(bootstrap)}$")
    plt.show()    
        
def booti(o_chain,N):
    boot_array=[]
    boot_array_final=[]
    for i in range(N):
        for j in range(len(o_chain)):
            boot_array.append(random.choice(o_chain))
        boot_array_final.append(np.average(boot_array))
        boot_array=[]
    return boot_array_final

def booti_ens(N):
    global N_cfg
    booti_error_ens=[]
    for N_cfg in np.arange(50000,1050000,50000):
        o_chain=markov_chain()
        o_chain_new=binning(o_chain,N)
        o_chain_boot=booti(o_chain_new, 200)
        booti_error_ens.append(np.sqrt(np.var(o_chain_boot)))
    plt.plot(np.arange(50000,1050000,50000),booti_error_ens)
    plt.xlabel(r"$N_{ens}$")
    plt.ylabel(r"$\sqrt{var(bootstrap)}$")
    plt.show()
        
        
        
#leapfrog_plot()
N_md=3
o_chain=markov_chain()
print(np.sqrt(np.var(o_chain)))
tau_int=int_auto_corr(o_chain)
# print(tau_int,'tau_int_before')
# print(np.average(o_chain))
show_auto_dec(o_chain,tau_int)
show_error_inc(o_chain,tau_int)
#print_nor_auto(o_chain,tau_int,str(1),str(2))
booti_ens(20)






        
    
