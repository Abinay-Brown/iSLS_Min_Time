'''
Optimization Parameters

'''
import numpy as np

N = 10
nx = 4
nu = 2
#tau = np.linspace(0, 1, N)
x0 = np.array([2, 2, 0, 0])            # initial conditions
Tf = 5
dt = Tf/(N-1)
xlim = 10
vlim = 1.5
ulim = 0.5
Q   = np.diag([1, 1, 1, 1])
R   = np.diag([1, 1]) 

#xf = np.array([7.2323, 5.4232, 0, 0])  # final conditions
#Exo = np.diag([0.1, 0.1, 0.005, 0.005]) # Exogenous Disturbance
