# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""

import numpy as np

L = 100 #meter
T = 
N = 100 # N + 1 elementer
dt = 
dN = L / N

# 
def K(z):
  K_0 = 1
  K_a = 1
  K_b = 1
  z_a = 1
  z_b = 1
  return K_0 + K_a*(z/z_a)*np.exp(-z/z_a) + (K_b*(L-z)/z_b) * np.exp(-(-L-z)/z_b)

# Konstrueerer diagonal matrise

L = diags([L_upper, L_main, L_lower], offsets = [1, 0, -1])
R = diags([R_upper, R_main, R_lower], offsets = [1, 0, -1])
