import numpy as np
import casadi as cas
from Parameters import *

# Dynamics
Ak = np.array([[1, 0, dt, 0], [0, 1, 0, dt], [0, 0, 1, 0], [0, 0, 0, 1]]).reshape((nx, nx))
Bk = np.array([[dt*dt/2, 0], [0, dt*dt/2], [dt, 0], [0, dt]])
Ik = np.eye(nx)

class fast_sls:
    def __init__(self):
        return True
    
    def init_dec_var(self):
        # Decision Variables
        self.z = cas.MX.sym('z', nx, N)
        self.v = cas.MX.sym('v', nu, N-1)
        return True

    def update_tightening(self, eps):
        # Tightening Parameters
        self.beta = np.zeros((N, N))
        self.hctk = np.zeros((N, 1))
        for i in range(0, N):
            for j in range(0, i+1):
                self.hctk[i] = np.sqrt(beta[i, j] + eps)
        return True
    
    def update_betas(self):
        
        return True
    def init_constraints(self):
        # Construct Equality Constraints
        Aeq = np.zeros(((N-1)*nx + nx, N*nx + (N-1)*nu))
        Aeq[0:nx, 0:nx] = np.eye(nx)
        offset = 0
        lbeq = np.zeros(((N-1)*nx + nx, N*nx + (N-1)*nu))
        ubeq = np.zeros(((N-1)*nx + nx, N*nx + (N-1)*nu))

        for i in range(N-1):
            r0 = nx + i*nx
            r1 = r0 + nx
            Aeq[r0 : r1, offset : (offset + nx)] = Ak
            Aeq[r0 : r1, (offset + nx)  : (offset + nx + nu)] = Bk
            Aeq[r0 : r1, (offset + nx + nu)  : (offset + nx + nu + nx)] = -Ik
            offset = offset + (nx + nu)
# Tightening Parameters
beta = np.zeros((N, N))
for i in range()

# Construct Inequality Linear Constraints    
    