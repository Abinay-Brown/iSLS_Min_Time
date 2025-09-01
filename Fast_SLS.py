import numpy as np
import casadi as cas
from Parameters import *

# Dynamics
Ak = np.array([[1, 0, dt, 0], [0, 1, 0, dt], [0, 0, 1, 0], [0, 0, 0, 1]]).reshape((nx, nx))
Bk = np.array([[dt*dt/2, 0], [0, dt*dt/2], [dt, 0], [0, dt]])
Ik = np.eye(nx)

class fast_sls:
    def __init__(self):
        return
    
    def init_dec_var(self):
        # Decision Variables
        self.z = cas.MX.sym('z', nx, N)
        self.v = cas.MX.sym('v', nu, N-1)
        return True
    def init_parameters(self, eps = 0):
        # Each row of Gk is gk,1...nc, selects each state and control variables
        nz = nx + nu  # Combined state and control dimensions
        nc = 2*nz     # Num of Constraints at each time step (upper and lower bound)
        self.Gk = np.zeros((nc, nz, N))
        for i in range(N):
            self.Gk[0  : nz, 0:nz, i] = - np.eye(nz)
            self.Gk[nz : nc, 0:nz, i] =   np.eye(nz)
        
        self.Phi_x = np.zeros((N * nx, N * nx))
        self.Phi_u = np.zeros((N * nu, N * nx))
        
        # Initialize Tightenings
        # Loop through row blocks
        for k in range(N):
            # Loop through column blocks
            hct_k = np.zeros((nc, 1))
            beta_k = np.zeros((nc, k+1))
            for j in range(k+1):
                Phi_kj = np.vstack((self.Phi_x[k*nx: (k+1)*nx, j*nx: (j+1)*nx], self.Phi_u[k*nu: (k+1)*nu, j*nx: (j+1)*nx]))
                # Loop through each constraint
                H_kj = np.zeros((nc, 1))
                for i in range(nc):
                    gkT = self.Gk[i, :, k]
                    H_kj[i, 0] = np.linalg.norm(gkT@Phi_kj)**2
                hct_k = hct_k + np.sqrt(H_kj + eps)
            print(beta_k.shape)
            beta_k[:, k] = hct_k.reshape((nc)) 


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

optim = fast_sls()
optim.init_dec_var()
optim.init_parameters()