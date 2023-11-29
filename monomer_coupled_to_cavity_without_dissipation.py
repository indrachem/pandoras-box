# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 18:28:39 2023

@author: Indranil
Coupling without dissipation (Monomer)

Parameters taken from Nat Comm (2023) 14:4804; CKT Weatherly, J Provazza, E A Weiss, R Tempelaar
"""

import numpy as np
import matplotlib.pyplot as plt

omega_c = 1600 #cm-1
omega_v = 1600 #cm-1
E_exc = 0 #cm-1 (Not concerned abt the energy-gap from ground state)
hbar = 6.626e-34/(2*np.pi)

s = 0.75 #Huang-Rhys factor
d = (2*s)*0.5 #s=0.5*(d)^2

g_val = np.linspace(0,0.9*219.474,5) #1 mHa = 219.474 cm-1

nvib = 8
order = 2*nvib**2
hf = int(order/2)
theta = 0.1 #Kelvin

from itertools import product
lst = list(product([j for j in range(nvib)], repeat=2)) #repeat = num of oscillators
lst = [list(i) for i in lst]

cav_plot = []
time_plot = []

for ctr in range(len(g_val)):

    H = np.zeros([order,order],float)
    
    for i in range(hf): # For diagonal elements (mol+cav)
        # index1 = int(i/nvib) #For nvib = 2, index1 = 0, 0, 1, 1
        # index2 = i - (index1*nvib) #index2 = 0, 1, 0, 1
        index1 = lst[i][0]
        H[i][i] = E_exc + omega_v*(index1+0.5)# E_vib = (n+0.5)*h_bar*w
        H[i+hf][i+hf] = omega_c*(index1+0.5)
    
    for i in range(hf): # Linear Vibronic Coupling Part
        # index1r = int(i/nvib) #For nvib = 2, index1r = 0, 0, 1, 1
        # index2r = i - (index1r*nvib) #index2r = 0, 1, 0, 1
        index1r = lst[i][0] 
        index2r = lst[i][1]
        for j in range(hf):
            index1c = lst[j][0]
            index2c = lst[j][1]
            # index1c = int(j/nvib) #index1c = 0, 0, 1, 1
            # index2c = j-(index1c*nvib) #index2c = 0, 1, 0, 1
            if index2r == index2c and index1r-index1c == 1: ##mol domain
                v = np.sqrt((index1c+1)/2)
                H[i][j] = -(omega_v)*d*v
                H[j][i] = H[i][j]
    
    # Ignoring Dipole Self Energy, Polarized Fock State Terms
    # H_int ~ H_cav-vib = (2g/h_bar)*((omega_c*omega_v)**0.5)*Q*q_c
    
    for i in range(hf,2*hf):
        # index1r = int(i/nvib) #For nvib = 2, index1r = 0, 0, 1, 1
        # index2r = i - (index1r*nvib) #index2r = 0, 1, 0, 1
        index1r = lst[i-hf][0] 
        index2r = lst[i-hf][1]
        for j in range(0,hf):
            index1c = lst[j][0]
            index2c = lst[j][1]
            # index1c = int(j/nvib) #index1c = 0, 0, 1, 1
            # index2c = j-(index1c*nvib) #index2c = 0, 1, 0, 1
            if index1r - index1c == 1 and index2r-index2c == 1:
                v = np.sqrt((index1c+1)/2)
                w = np.sqrt((index2c+1)/2)
                H[i][j] = (2*g_val[ctr])*v*w
                H[j][i] = H[i][j]
                
                H[i-hf][j+hf] = (2*g_val[ctr])*v*w
                H[j+hf][i-hf] = H[i-hf][j+hf]
    
    #%%     
    eigVals, eigVecs = np.linalg.eigh(H)
    idx = eigVals.argsort()[::1] #gives the index of the elements to be sorted
    eigVals = eigVals[idx] #sorts acc to the index idx
    eigVecs = eigVecs[:,idx]
    
    # Population Dynamics
    
    time = []
    mol = []
    cav = []
    P_mol = []
    P_cav = []
    
    xp = 1
    yp = 0
    c = 2*np.pi*2.99792458*(10**10)
    
    s = 0
    for t in np.arange(0,400*(10**(-15)),5*(10**(-15))):
        mol.append([])
        cav.append([])
    
        for i in range(0,int(order/2)):
            s1 = 0
            s2 = 0
    
            for k in range(0,order):
                s1 = s1 + (eigVecs[0][k]* xp + eigVecs[int(order/2)][k]* yp)*(np.exp(-1j *t*c*eigVals[k])) *eigVecs[i][k]
                s2 = s2 + (eigVecs[0][k]* xp + eigVecs[int(order/2)][k]* yp)*(np.exp(-1j *t*c*eigVals[k])) *eigVecs[i+int(order/2)][k]
           
            mol[s].append((abs(s1)**2))
            cav[s].append((abs(s2)**2))
            P1 = sum(mol[s])
            P2 = sum(cav[s])
        time.append(t)
       
    
        P_mol.append(P1)
        P_cav.append(P2)
        s+= 1
    
    P_mol = np.array(P_mol)
    P_cav = np.array(P_cav)
    
    time_plot.append(time)
    cav_plot.append(P_cav)
    
    # Transfer=100*(np.max(P_mol)/(P_mol[0]+P_cav[0]))
    
#%%
fig, axes = plt.subplots(len(g_val), sharex=True, sharey=True, figsize=(8,8))
for i in range(len(g_val)):
    axes[i].plot(time_plot[i],cav_plot[i])
    if i == len(g_val)-1:
        axes[i].set_xlabel('time / ps')
    axes[i].set_ylabel('population')
    axes[i].set_title(f'g = {g_val[i]} cm-1')

plt.xticks([0,1.0e-13,2.0e-13,3.0e-13,4.0e-13],[0,0.1,0.2,0.3,0.4])
# plt.yticks([0,0.5,1])

plt.tight_layout()

plt.show()