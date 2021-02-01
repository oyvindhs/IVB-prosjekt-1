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

# 
def K(z):
  K_0 = 1
  K_a = 1
  K_b = 1
  z_a = 1
  z_b = 1
  return K_0 + K_a*(z/z_a)*np.exp(-z/z_a) + (K_b*(L-z)/z_b) * np.exp(-(-L-z)/z_b)

# Konstrueerer diagonal matrise

alpha = dt / (2 * dz**2)
gamma = 2*alpha*k_w*dz*(1 - (1 - (K_1 - K_0) / (2*K_0))
K_array = np.zeros(N + 1)


L_upper = np.zeros(N-1)
L_main = np.zeros(N)
L_lower = np.zeros(N-1)

L = diags([L_upper, L_main, L_lower], offsets = [1, 0, -1])
R = diags([R_upper, R_main, R_lower], offsets = [1, 0, -1])

 print("hallo")



