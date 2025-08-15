'''
Optimization Parameters

'''
import numpy as np

N = 50
nx = 4
nu = 2
tau = np.linspace(0, 1, N)
x0 = np.array([0, 0, 0, 0])            # initial conditions
xf = np.array([7.2323, 5.4232, 0, 0])  # final conditions
xlim = 10
vlim = 1
ulim = 0.5
Exo = np.diag([0.1, 0.1, 0.05, 0.05]) # Exogenous Disturbance

