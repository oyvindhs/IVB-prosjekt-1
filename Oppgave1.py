# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""

import numpy as np

L = 100 #meter
T = 1
N = 100 # N + 1 elementer
dt = 1
dz = L / N

# Diffusjonsparameter
def K(z):
    
    K_0 = 10**(-3) #m^2
    K_a = 2*10**(−2) #m^2
    z_a = 7 #m
    K_b = 5 * 10**(−2) #m^2
    z_b = 10 #m
  
  return K_0 + K_a*(z/z_a)*np.exp(-z/z_a) + (K_b*(L-z)/z_b) * np.exp(-(-L-z)/z_b)

# Konstruerer array med diffusjonsjonsparameter K_n for vann-element n
K_ray = K(np.linspace(0, L, N+1))



# Konstrueerer tridiagonal matrise

alpha = dt / (2 * dz**2)
gamma = 2*alpha*k_w*dz*(1 - (1 - (K_1 - K_0) / (2*K_0)))
K_merket = np.zeros(N-2)
for m in range(1, N): #OBS! Riktig lengde på K_merket?
    K_merket[m] = K_ray[m+1] - K_ray[m-1] 
  

L_upper = np.zeros(N-1)


L_main = np.zeros(N)

L_lower = np.zeros(N-1)

L = diags([L_upper, L_main, L_lower], offsets = [1, 0, -1])



