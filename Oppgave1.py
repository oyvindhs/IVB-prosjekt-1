# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""

import numpy as np
from scipy.sparse import diags
from matplotlib import pyplot as plt

L = 100 #meter
T = 180 #år
N = 100 # N + 1 elementer
dt = 1
dz = L / N
k_w = 6.97*10**(-5) #m/s
c_eq = 2.1 #mol/m^3

# Diffusjonsparameter
def K(z):
    
    K_0 = 10**(-3) 
    K_a = 2*10**(-2) 
    z_a = 7 #m
    K_b = 5 * 10**(-2) #m^2
    z_b = 10 #m
  
    return K_0 + K_a*(z/z_a)*np.exp(-z/z_a) + (K_b*(L-z)/z_b) * np.exp(-(-L-z)/z_b)

# Konstruerer array med diffusjonsjonsparameter K_n for vann-element n
K_ray = K(np.linspace(0, L, N+1))



# Konstrueerer tridiagonal matrise


alpha = dt / (2 * dz**2)
gamma = 2*alpha*k_w*dz*(1 - (K_ray[1] - K_ray[0]) / (2*K_ray[0]))
K_merket = np.zeros(N-2)
for m in range(1, N-2): #OBS! Riktig lengde på K_merket?
    K_merket[m] = K_ray[m+1] - K_ray[m-1] 
  

#Konstruerer de tre diagonalene i L-matrisen (av hhv. lengde N, N + 1, N), ut fra oppgitt matrise
L_upper = [-2*alpha*K_ray[0]]
for i in range(1, N):
    L_upper.append((-alpha/4)*K_merket[i] - alpha*K_ray[i])
  

L_mid = [1 + 2*alpha*K_ray[0] + gamma]
for i in range(1, N + 1):
    L_mid.append((1 + 2*alpha*K_ray[i]))
  

L_lower = [(alpha/4)*K_merket[1] - alpha*K_ray[1]]
for i in range(2, N):
    L_lower.append((alpha/4)*K_merket[i] - alpha*K_ray[i])
L_lower.append(-2*alpha*K_ray[N])

#Konstruerer de tre diagonalene i R-matrisen (av hhv. lengde N, N + 1, N), ut fra oppgitt matrise


R_upper = [+2*alpha*K_ray[0]]
for i in range(1, N):
    R_upper.append((+alpha/4)*K_merket[i] + alpha*K_ray[i])
 

R_mid = [1 - 2*alpha*K_ray[0] - gamma]
for i in range(1, N + 1):
    R_mid.append((1 - 2*alpha*K_ray[i]))
  

R_lower = [(-alpha/4)*K_merket[1] + alpha*K_rat[1]]
for i in range(2, N):
    R_lower.append((-alpha/4)*K_merket[i] + alpha*K_rat[i])
R_lower.append(2*alpha*K_ray[N])


L = diags([L_upper, L_mid, L_lower], offsets = [1, 0, -1])
R = diags([R_upper, R_mid, R_lower], offsets = [1, 0, -1])






def simulate():
    S_vec = np.zeros(N+1)
    S[0] = 2*gamma*c_eq




