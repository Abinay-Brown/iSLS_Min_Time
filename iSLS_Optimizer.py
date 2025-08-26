import scipy.io
import numpy as np
import cvxpy as cp
from Parameters import *


class SLS_optimizer:

    def __init__(self):

        return

    def initParameters(self):
        self.z = cp.Parameter((nx, N))
        self.v = cp.Parameter((nu, N))

        self.Ablk = cp.Parameter((N * nx, N * nx))  # A Block Matrix
        self.Bblk = cp.Parameter((N * nx, N * nu))  # B Block Matrix 
        self.Zblk = cp.Parameter((N * nx, N * nx))  # Z Block Matrix
        self.Iblk = cp.Parameter((N * nx, N * nx))  # I Block Matrix
        self.Qblk = cp.Parameter((N * nx, N * nx))  # State Weight Block Matrix 
        self.Rblk = cp.Parameter((N * nu, N * nu))  # Control Weight Block Matrix

        self.E    = cp.Parameter((nx, nx))          # Exogenous Disturbance Diagonal Matrix
        #self.mu   = cp.Parameter((nx, nx))          # Linearization Error Diagonal Matrix
        
        self.cx   = cp.Parameter((nx, nx))
        self.cu   = cp.Parameter((nu, nu))
        
        self.bx   = cp.Parameter((nx, 2))
        self.bu   = cp.Parameter((nu, 2))
        #self.tau   = cp.Parameter((N)) # Double check the dimension if N + 1
         
        return True

    def setParameters(self, Xnom, Unom, Tnom, mu):
        
        # Set Nominal Trajectory
        self.z.value = Xnom
        self.v.value = Unom
        
        # Set Block Matrices
        Amat, Bmat, Zmat = self.Linearize(Xnom, Unom, Tnom)
        self.Ablk.value = Amat
        self.Bblk.value = Bmat
        self.Zblk.value = Zmat
        self.Iblk.value = np.eye((N * nx))
        self.Qblk.value = np.kron(np.eye(N), Q)
        self.Rblk.value = np.kron(np.eye(N), R)

        # Set Disturbance Matrices
        self.E.value = Exo
        #self.mu.value = mu

        # Set c select
        self.cx.value = np.eye(nx)
        self.cu.value = np.eye(nu)
        
        # Set b for bounds
        self.bx.value = np.array([[xlim, -xlim], [xlim, -xlim], [vlim, -vlim], [vlim, -vlim]])
        self.bu.value = np.array([[ulim, -ulim], [ulim, -ulim]])
        #self.tau.value = np.zeros((N))
        return True
    def initDecisionVariables(self):
        #self.z = cp.Variable((nx, N))
        #self.v = cp.Variable((nu, N))
        self.Phi_x = cp.Variable((N * nx, N * nx))
        self.Phi_u = cp.Variable((N * nu, N * nx))
         
        return True
    def initConstraints(self):
        self.Constraints = []
        # Boundary Conditions
        
        # Phi_x & Phi_u upper diagonal constraint
        for i_blk in range(N):        
            for j_blk in range(N):    
                if j_blk > i_blk:     
                    self.Constraints += [self.Phi_x[i_blk*nx:(i_blk+1)*nx, j_blk*nx:(j_blk+1)*nx] == 0]
                    self.Constraints += [self.Phi_u[i_blk*nu:(i_blk+1)*nu, j_blk*nx:(j_blk+1)*nx] == 0]
    
                    
        # [I - ZA - ZB]*Phi = [I] Constraint
        #print(P)
        #Phi = cp.vstack([self.Phi_x, self.Phi_u])
        #print(Phi.shape)
        #print(((self.Iblk - self.Zblk@self.Ablk - self.Zblk@self.Bblk)@Phi).shape)
        
        self.Constraints += [((self.Iblk - self.Zblk@self.Ablk)@self.Phi_x - (self.Zblk@self.Bblk@self.Phi_u)) == self.Iblk]
        
        # State Tube Constraint
        # Loop through the Trajectory points
        for i in range(N):
            # Loop through each state
            for j in range(nx):
                # Loop to index block matrices
                LHS_upper = 0
                LHS_lower = 0
                for k in range(i+1):
                    LHS_upper = LHS_upper + cp.norm(self.cx[j, :]@self.Phi_x[i*nx: (i+1)*nx, k*nx: (k+1)*nx]@self.E, 'inf')
                    #LHS_lower = LHS_lower + cp.norm(-self.cx[j, :]@self.Phi_x[i*nx: (i+1)*nx, k*nx: (k+1)*nx]@self.E, 'inf')               
                # Set Upper Tube Bound
                LHS_upper = LHS_upper + self.cx[j, :]@self.z[:, i] - self.bx[j, 0]     
                # Set Lower Tube Bound  
                LHS_lower = LHS_lower - self.cx[j, :]@self.z[:, i] + self.bx[j, 1] # Check if it is -bx !!!!!!!!!!!!!!!!!
                self.Constraints += [LHS_upper <= 0]
                self.Constraints += [LHS_lower <= 0]
        
        # Control Tube Contraint
        # Loop through the Trajectory points
        for i in range(N):
            # Loop through each state
            for j in range(nu):
                # Loop to index block matrices
                LHS_upper = 0
                LHS_lower = 0
                for k in range(i+1):
                    LHS_upper = LHS_upper + cp.norm(self.cu[j, :]@self.Phi_u[i*nu: (i+1)*nu, k*nx: (k+1)*nx]@self.E, 'inf')
                    LHS_lower = LHS_lower + cp.norm(-self.cu[j, :]@self.Phi_u[i*nu: (i+1)*nu, k*nx: (k+1)*nx]@self.E, 'inf')               
                # Set Upper Tube Bound
                LHS_upper = LHS_upper + self.cu[j, :]@self.v[:, i] - self.bu[j, 0]
                # Set Lower Tube Bound  
                LHS_lower = LHS_lower - self.cu[j, :]@self.v[:, i] + self.bu[j, 1] # Check if it is -bu !!!!!!!!!!!!!!!!!
                self.Constraints += [LHS_upper <= 0]
                self.Constraints += [LHS_lower <= 0]
        '''
        # 27e constraint
        
        for i in range(1, N):
            LHS = 0
            for k in range(i):
                LHS = LHS + cp.norm(self.Phi_x[(i-1)*nx: i*nx, k*nx: (k+1)*nx]@cp.hstack([self.E, self.mu*self.tau[k]**2]), 'inf')               
            self.Constraints += [LHS <= self.tau[i]]
        
        for i in range(1, N):
            LHS = 0
            for k in range(i):
                LHS = LHS + cp.norm(self.Phi_u[(i-1)*nu: i*nu, k*nx: (k+1)*nx]@cp.hstack([self.E, self.mu*self.tau[k]**2]), 'inf')               
            self.Constraints += [LHS <= self.tau[i]]
        '''    
        return True

    def solve_iSLS(self):
        cost = 0

        for i in range(1, N):
            sub = 0
            for k in range(i):
                sub = sub + self.Phi_x[(i-1)*nx: i*nx, k*nx: (k+1)*nx]@self.E               
            cost = cost + cp.norm(self.Qblk[(i-1)*nx: i*nx, (i-1)*nx: i*nx]@sub, 'fro')
        
        for i in range(1, N):
            sub = 0
            for k in range(i):
                sub = sub + self.Phi_u[(i-1)*nu: i*nu, k*nx: (k+1)*nx]@self.E               
            cost = cost + cp.norm(self.Rblk[(i-1)*nu: i*nu, (i-1)*nu: i*nu]@sub, 'fro')
        
        self.Phi_x.value = np.zeros((N * nx, N * nx))
        self.Phi_u.value = np.zeros((N * nu, N * nx))
        prob = cp.Problem(cp.Minimize(cost), self.Constraints)
        prob.solve(solver=cp.ECOS, verbose=True,
           max_iters=5000, warm_start=True)
        return True
    def computeLinearizationError(self):
        mu = 0
        return mu

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
Xnom = Nominal['Xnom']
Unom = Nominal['Unom']
Tnom = Nominal['Tnom'].item()    
Xnom = Xnom[:, 0:-1]

opt = SLS_optimizer()
opt.initParameters()
opt.setParameters(Xnom, Unom, Tnom, np.diag([0, 0, 0, 0]))
opt.initDecisionVariables()
opt.initConstraints()
opt.solve_iSLS()
 