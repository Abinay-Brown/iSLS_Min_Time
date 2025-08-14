import scipy.io
import numpy as np
import cvxpy as cp
from Parameters import *


class SLS_optimizer:

    def __init__(self):

        return

    def initParameters(self, Xnom, Unom):
        Amat, Bmat, Zmat = self.Linearize(Xnom, Unom)
        self.Ablk = cp.Parameter((N * nx, N * nx))
        self.Bblk = cp.Parameter((N * nx, N * nu))
        self.Zblk = cp.Parameter((N * nx, N * nx))
        self.Iblk = cp.Parameter((N * nx, N * nx))

        self.Ablk.value = Amat
        self.Bblk.value = Bmat
        self.Zblk.value = Zmat
        self.Iblk.value = np.eye((N * nx))

        # Need to setup E matrix
        # Need to setup mu: for that solve hessian maximization problem
        # Need to setup c to select elements
        # Need to setup b for bounds
        return True

    def initDecisionVariables(self):
        self.z = cp.Variable((nx, N))
        self.v = cp.Variable((nu, N))
        
        self.Phi_x = cp.Variable((N * nx, N * nx))
        self.Phi_u = cp.Variable((N * nu, N * nx))
        self.tau   = cp.Variable((N)) # Double check the dimension if N + 1
         
        return True
    def initConstraints(self):
        self.Constraints = []
        # Boundary Conditions
        # Dynamics Constraints
        # Control limits
        
        # Phi_x & Phi_u upper diagonal constraint
        for i_blk in range(N):        
            for j_blk in range(N):    
                if j_blk > i_blk:     
                    self.constraints += [Phi_x[i_blk*nx:(i_blk+1)*nx, j_blk*nx:(j_blk+1)*nx] == 0]
                    self.constraints += [Phi_u[i_blk*nu:(i_blk+1)*nu, j_blk*nx:(j_blk+1)*nx] == 0]
    
                    
        # [I - ZA - ZB]*Phi = [I] Constraint
        Phi = cp.vstack([self.Phi_x, self.Phi_u])
        self.Constraints += [(self.Iblk - self.Zblk@self.Ablk - self.Zblk@self.Bblk)@phi == self.Iblk]
        # 27d constraint 
        # 27e constraint
            
        return
    
    def updateParameters(self):
        
        return True 

    def Linearize(self, Xnom, Unom, Tnom):
        Ablk = np.zeros((N * nx, N * nx))
        Bblk = np.zeros((N * nx, N * nu))
        Zblk = np.zeros((N * nx, N * nx))
        idx = nx 
        idu = nu
        dt = Tnom/N
        for i in range(N):
            Ablk[idx-nx : idx, idx-nx : idx] = np.array([[1, 0, dt, 0], [0, 1, 0, dt], [0, 0, 1, 0], [0, 0, 0, 1]]).reshape((nx, nx))
            Bblk[idx-nx : idx, idu-nu : idu] = np.array([[dt*dt/2, 0], [0, dt*dt/2], [dt, 0], [0, dt]])
            if i > 0:
                Zblk[idx-nx:idx, idx-2*nx:idx-nx] = np.eye(nx)
            idx = idx + nx
            idu = idu + nu
        return Ablk, Bblk, Zblk
    
Nominal = scipy.io.loadmat('Guess_Trajectory.mat')
Xnom = Nominal['Xnom'].T
Unom = Nominal['Unom'].T
Tnom = Nominal['Tnom'].item()    
opt = SLS_optimizer()
Ablk, Bblk, Zblk = opt.Linearize(Xnom, Unom, Tnom)
 