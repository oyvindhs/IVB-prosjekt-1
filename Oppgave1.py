
import numpy as np
from scipy.sparse import diags
from matplotlib import pyplot as plt
from scipy.integrate import simps
from numba import jit

CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
color_list = [CB91_Blue, CB91_Pink, CB91_Green, CB91_Amber,
              CB91_Purple, CB91_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)


L = 100 #meter
T = 180 * 24 * 60 * 60  #sekunder
N = 1000 # N + 1 elementer
dt = 60*60
dz = L / N
a = 6.97e-7 #s/m
u = 10 #m/s
k_w = a*u**2  #m/s, (6.97e-5)
p_co2 = 415e-6 #atm
H = 5060 #mol/(m^3atm)
c_eq = H*p_co2 #mol/m^3 (ca. 2.1)

# Diffusjonsparameter
def K(z):
    
    K_0 = 10**(-3) 
    K_a = 2*10**(-2) 
    z_a = 7 #m
    K_b = 5 * 10**(-2) #m^2
    z_b = 10 #m
  
    return K_0 + K_a*(z/z_a)*np.exp(-z/z_a) + (K_b*(L-z)/z_b) * np.exp(-(L-z)/z_b)

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

def iterate(C_vec, R_matrix, S_vec, L_matrix):
    V_matrix = R_matrix.dot(C_vec) + S_vec
    #iterate C by solving L C_next = V
    C_vec = tdma(L_matrix, V_matrix)
    
    return C_vec

def simulate():
    C_timeline = []
    S_vec = np.zeros(N+1)
    S_vec[0] = 2*gamma*c_eq
    C_vec = np.zeros(N+1)
    
    for t in range(int(T/dt)):
        if (t%(10**(7)) == 0):
            print(round(t/T, 3)*100, "%")
        C_timeline.append(C_vec)
        C_vec = iterate(C_vec, R_matrix, S_vec, L_matrix)    
    print("100%")
    return C_timeline # [tidspunkt][koordinat]

    
output = simulate()  

t = np.linspace(0, T, len(output)) 

max_ar = np.zeros(len(output))
min_ar = np.zeros(len(output))
for i in range(1, len(output)):
    max_ar[i] = max(output[i])
    min_ar[i] = min(output[i])


plt.figure(figsize=(12,8))
plt.plot(t, max_ar)
plt.plot(t, min_ar)
plt.show()

plt.figure(figsize=(12,6))
#for i in range(0, len(output),):
plt.plot(output[0], -np.linspace(0, L, N + 1)) #t = 0
plt.plot(output[1], -np.linspace(0, L, N + 1)) #t = 60*60s = 1h
plt.plot(output[2], -np.linspace(0, L, N + 1))
plt.plot(output[3], -np.linspace(0, L, N + 1))
plt.plot(output[4], -np.linspace(0, L, N + 1))
plt.plot(output[5], -np.linspace(0, L, N + 1))
plt.plot(output[6], -np.linspace(0, L, N + 1))
plt.plot(output[7], -np.linspace(0, L, N + 1))
plt.plot(output[10], -np.linspace(0, L, N + 1)) #t = 60*60*10s = 10h
plt.plot(output[50], -np.linspace(0, L, N + 1)) #t = 50h
plt.plot(output[500], -np.linspace(0, L, N + 1)) #t = 500h = ca 21d
plt.plot(output[2000], -np.linspace(0, L, N + 1)) #t = 2000h = ca 83d
plt.plot(output[len(output)-1], -np.linspace(0, L, N + 1)) #t = 180d
plt.show()
