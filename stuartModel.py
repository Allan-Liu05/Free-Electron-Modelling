import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
# Constants
pulse_duration = 150e-15
pulse_shift = 0
I_peak = 60 # Peak intensity TW/cm^2
P_constant = 9.52e10  # Constant in P(I) 
tScale = 175 # half the timescale of which the model is to simulate

# Keldysh Parameter for out laser is 1.05448077238 (100% in multiphoton regime)


# Post Laser Pulse Constants
#alpha = 8.597e-12
alpha = 1.4e-8
beta = 2.5940830498e-31

# Define the Gaussian function for I(t)
def I_t(t):
    return I_peak*np.exp(-4*math.log(2)*(((t-pulse_shift)/pulse_duration)**2))

def I_t2(t):
    return 12*np.exp(-4*math.log(2)*(((t-200e-15)/100e-15)**2))

def I_t3(t):
    return 3.5*np.exp(-4*math.log(2)*(((t-1000e-15)/1000e-15)**2))

# Define P(I) as a function of time
def P_I(t):
    I = I_t(t)
    return P_constant * I**8 * 1e12 
    
def avalanche(t):
    I = I_t(t)
    return I*1e12


def rate_equation(n,t):

    P = P_I(t)
    ava = avalanche(t)
    ava_efficiency = 1 - n  /(6.6e22*(8/3))
    mpi_efficiency = 1
    A = 1 - n/6.6e22

    if n >= 6.6e22*(8/3):
        return 6.6e22*(8/3) 

    if n > 6.6e22:
        mpi_efficiency = .25
    else: 
        mpi_efficiency = 1

    print(A*P)
    return mpi_efficiency*A*P + ava_efficiency*4*n*ava - n/(150e-15)

def rateTest(n,t):
    
    P = P_I(t)
    ava = avalanche(t)
    ava_efficiency = 1 - n  /(6.6e22*(8/3))
    mpi_efficiency = 1
    A = 1 - n/6.6e22

    if n >= 6.6e22*(8/3):
        return 6.6e22*(8/3) 

    if n > 6.6e22:
        mpi_efficiency = .25
    else: 
        mpi_efficiency = 1

    print(A*P)
    return mpi_efficiency*A*P -n/(90e-15)#+ ava_efficiency*4*n*ava - n/(90e-15)

N_e = 1.327e23
neutral_atoms = 0
def postLaserRate(n, t):
    global N_e
    global neutral_atoms
    positive_ions = (3/8)*n
    neutral_atoms += (3/8)*abs(n-N_e)
    N_e = n
    return alpha*neutral_atoms*n - beta*positive_ions*(n**2)


t = np.linspace(-tScale,tScale,1000)*1e-15

# plt.plot(t, I_t(t))
# plt.show()

# n=0
# y = []
# for dt in t:
#     n += dt*(P_I(dt)+n*avalanche(dt))
#     y.append(n)
# plt.plot(t,y)
# plt.yscale('log')
# plt.show()

n_t = odeint(rate_equation,0, t)
n_t2 = odeint(rateTest, 0, t)

fig, ax1 = plt.subplots()
ax1.plot(t*1e15, n_t, color='b', label="Multiphoton and Avalanche Ionization")
ax1.plot(t*1e15, n_t2, color='r', label="Multiphon Ionization Only")
ax1.set_xlabel("Time (fs)")
ax1.set_ylabel("Free Electron Density (cm^-3)")
ax1.set_yscale("log")
plt.legend()
# ax2 = ax1.twinx()
# ax2.plot(t*1e15, I_t(t), color='g', label="Gaussian Pulse")
# ax2.set_ylabel("Intensity (TW/cm^2)")
plt.show()

# n_t = solve_ivp(rate_equation,[min(t),max(t)], y0=[0],method="RK45", t_eval=t)
# sol = n_t.y[0]
# plt.plot(t*1e15,sol,label="Multiphoton and Avalanche Ionization")
# plt.yscale('log')
# plt.legend()
# plt.show()

# us = []
# s12 = []
# s35 = []
# us_add = 0
# s12_add = 0
# s35_add = 0
# t = np.linspace(0,2000,10000)*1e-15
# for dt in t:
#     us_add += (dt*I_t(dt))
#     s12_add += (dt*I_t2(dt))
#     s35_add += (dt*I_t3(dt))
#     us.append(us_add)
#     s12.append(s12_add)
#     s35.append(s35_add)

# plt.plot(t, us, color='r',label='Our Laser 100TW/cm^2')
# plt.plot(t, s12, color='g',label='Stuart 12TW/cm^2')
# plt.plot(t, s35,color='b', label="Stuart 3.5 TW/cm^2")
# plt.legend()
# plt.show()

# POST LASER STUFF
t_post = np.linspace(300,5000,1000)*1e-15
post_n_t = odeint(postLaserRate, 1.327e23, t_post)
plt.plot(t_post*1e15,post_n_t)
plt.xlabel("Time (fs)")
plt.ylabel("Free electron density (cm^-3)")
plt.yscale("log")
plt.show()