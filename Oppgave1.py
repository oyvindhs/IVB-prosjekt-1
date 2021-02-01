# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

L = 1

def K(z):
  K_0 = 1
  K_a = 1
  K_b = 1
  z_a = 0
  return K_0 + K_a*(z/z_a)*np.exp(-z/z_a) + (K_b*(L_z)/z_b) * np.exp(-(-L-z)/z_b)
