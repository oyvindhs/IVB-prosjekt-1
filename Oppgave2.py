# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:52:07 2021

@author: oy019
"""


import numpy as np
from scipy.sparse import diags
from matplotlib import pyplot as plt
from scipy.integrate import simps
from numba import jit


L = 4000 #meter
T = 10*365*24*60*60 #sek
N = 40000 # N + 1 elementer
dt = 60*60*24
dz = L / N
k_w = 6.97*10**(-5) #m/s


#Funksjon for av co2 i atmosfæren. Returnerer C_eq(tid)
def c_eq(t):
    
    return (415 + 2.3 * (t/(365*24*60*60))) * 5060*10**(-6)
    

# Diffusjonsparameter
def K(z):
    
    K_0 = 10**(-4) #m^2/s
    K_1 = 10**(-2) #m^2/s
    a = 0.5 #m^-1
    z_0 = 100 #m
  
    return K_1 + (K_0 - K_1) / (1 + np.exp(-a*(z - z_0)))

# Konstruerer array med diffusjonsjonsparameter K_n for vann-element n
K_ray = K(np.linspace(0, L, N+1))



# Konstrueerer tridiagonal matrise


alpha = dt / (2 * dz**2)
gamma = 2*alpha*k_w*dz*(1 - (K_ray[1] - K_ray[0]) / (2*K_ray[0]))
K_merket = np.zeros(N)
for m in range(1, N): #OBS! Riktig lengde på K_merket?
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

R_upper = [2*alpha*K_ray[0]]
for i in range(1, N):
    R_upper.append((alpha/4)*K_merket[i] + alpha*K_ray[i])
  
R_mid = [1 - 2*alpha*K_ray[0] - gamma]
for i in range(1, N + 1):
    R_mid.append((1 - 2*alpha*K_ray[i]))
  
R_lower = [(-alpha/4)*K_merket[1] + alpha*K_ray[1]]
for i in range(2, N):
    R_lower.append((-alpha/4)*K_merket[i] + alpha*K_ray[i])
R_lower.append(2*alpha*K_ray[N])



L_matrix = diags([L_upper, L_mid, L_lower], offsets = [1, 0, -1])
R_matrix = diags([R_upper, R_mid, R_lower], offsets = [1, 0, -1])

#Fra hjelpekode på BB:
#Løser matriseligning
@jit(nopython = True)
def tdma_solver(a, b, c, d):
    N = len(d)
    c_ = np.zeros(N-1)
    d_ = np.zeros(N)
    x  = np.zeros(N)
    c_[0] = c[0]/b[0]
    d_[0] = d[0]/b[0]
    for i in range(1, N-1):
        c_[i] = c[i]/(b[i] - a[i-1]*c_[i-1])
    for i in range(1, N):
        d_[i] = (d[i] - a[i-1]*d_[i-1])/(b[i] - a[i-1]*c_[i-1])
    x[-1] = d_[-1]
    for i in range(N-2, -1, -1):
        x[i] = d_[i] - c_[i]*x[i+1]
    return x

def tdma(A, b):
    # Solves Ax = b to find x
    x = tdma_solver(A.diagonal(-1), A.diagonal(0), A.diagonal(1), b)
    return x

def iterate_C(C_vec, R_matrix, S_cur, S_next, L_matrix):
    V_matrix = R_matrix.dot(C_vec) + 0.5 * (S_cur + S_next)
    #iterate C by solving L C_next = V
    C_vec = tdma(L_matrix, V_matrix)
    
    return C_vec

def simulate():
    C_timeline = []
    S_cur = np.zeros(N+1)
    S_cur[0] = 2 * gamma * c_eq(0)
    S_next = np.zeros(N+1)
    C_vec = np.ones(N+1)*2.1 #mol/m^3
    for t in range(0, T, dt):
        
        if (t%(10**(6)) == 0):
            print(round(t/T, 3)*100, "%")
            
        
        S_next[0] = 2*gamma*c_eq(t+dt)
        
        C_timeline.append(C_vec)
        C_vec = iterate_C(C_vec, R_matrix, S_cur, S_next, L_matrix)
        
        S_cur = S_next
        
    print("\n Simulated ", T, " seconds in", T/dt, " intervals. \n")
    print("Depth: ", L, "consisting og ", N+1, "elements. \n")
    return C_timeline # [tidspunkt][koor"dinat]
#-------------------------------------------------------------------------------------------------------
#Task 1
"""Plot the concentration as a function of depth, for the times 0, 2.5, 5 and 10 years into the simulation.
Comment on the results, in light of the illustration in Fig. 2. Starting at eq. """

output = simulate()  

#OBS, magiske tall!
"""
time_zero = output[0]
time_twoandhalf = output[91] 
time_five = output[182]
time_ten = output[364]

plt.figure()
plt.plot(time_zero, -np.linspace(0, L, N+1))
plt.plot(time_twoandhalf, -np.linspace(0, L, N+1))
plt.plot(time_five, -np.linspace(0, L, N+1))
plt.plot(time_ten, -np.linspace(0, L, N+1))
plt.show()
"""
#-------------------------------------------------------------------------------------------------------
#Task 3
"""Plot the total mass of DIC in the global oceans, as a function of time, for the ye"ars 2020–2"030"""

mass_ray = np.zeros(len(output))
time_ray = np.linspace(0, T, len(output))
for i in range(len(output)):
    mass = simps(output[i],np.linspace(0, L, N+1) )
    mass_ray[i] = mass * 12 * 10**(-3) * 360 * 10**12 #kg per column x sea surface
"""
plt.figure()
plt.plot(time_ray, mass_ray)
plt.show()
"""
#-------------------------------------------------------------------------------------------------------
#Task 4

"""Find the amount of CO2 absorbed by the entire global ocean in a year by looking at the mass in the
water column at the start of the simulation, compared to the mass at the end of the simulation, and
take the average over the 10 years."""

val = (mass_ray[-1] - mass_ray[0])/10
print("Average increase in mass of DIC in one year: ", round(val/1000,0), " metric tons." )


