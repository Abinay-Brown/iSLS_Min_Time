import numpy as np
import casadi as cas
from Parameters import *
import matplotlib.pyplot as plt
import scipy.io

# Decision Variables
X = cas.MX.sym('X', nx, N+1)
U = cas.MX.sym('U', nu, N)
T = cas.MX.sym('sigma')


def f(x, u):
    px, py, vx, vy = x[0], x[1], x[2], x[3]
    ax, ay = u[0], u[1]
    return cas.vertcat(vx, vy, ax, ay)

# Runge-Kutta 4th Order Integrator
def rk4_step(x, u, h):
    k1 = f(x, u)
    k2 = f(x + (h/2) * k1, u)
    k3 = f(x + (h/2) * k2, u)
    k4 = f(x +  h    * k3, u)
    return x + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

# Shooting Function
def shoot(x, u, dt, n_sub=10):
    h = dt / n_sub
    x_temp = x
    for _ in range(n_sub):
        x_temp = rk4_step(x_temp, u, h)
    return x_temp

g_list = []
dt = T/N    
for k in range(N):
    xk = X[:, k]       
    uk = U[:, k]        
    xkp1 = shoot(xk, uk, dt)  
    x_next = X[:, k+1]  
    res_k = x_next - xkp1 
    g_list.append(res_k)

g_cont = cas.vertcat(*g_list)

g_bc = cas.vertcat(X[:, 0] - x0, X[:, N] - xf)

z   = cas.vertcat(cas.vec(X), cas.vec(U), T)
nz  = int(z.numel())
lbx = [-cas.inf]*nz
ubx = [ cas.inf]*nz

# Offsets
offX = 0
offU = (N+1)*nx
idxT = offU + N*nu

for k in range(N+1):
    lbx[offX + k*nx + 0] = 0.0   # px min
    ubx[offX + k*nx + 0] = 10.0  # px max
    lbx[offX + k*nx + 1] = 0.0   # py min
    ubx[offX + k*nx + 1] = 10.0  # py max
u_max = 2.0
for k in range(N):
    for j in range(nu):
        idx = offU + k*nu + j
        lbx[idx] = -u_max
        ubx[idx] =  u_max

T_min, T_max = 1e-2, 100.0
lbx[idxT], ubx[idxT] = T_min, T_max

J = T +  (T/N)*cas.sumsqr(U)

# --- Combine all constraints ---
g = cas.vertcat(g_bc, g_cont)
lbg = [0]*g.numel()
ubg = [0]*g.numel()

# --- Build and solve NLP ---
nlp = {'x': z, 'f': J, 'g': g}
opts = {"ipopt.print_level": 5, "ipopt.sb": "yes", "print_time": 1}  # verbose iterations
solver = cas.nlpsol('solver', 'ipopt', nlp, opts)

# Initial guess
z0 = [0]*(nz)
z0[idxT] = 2.0  # positive initial time guess (adjust as you like)

# Solve
sol = solver(x0=z0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
z_opt = sol['x'].full().flatten()

''' Extract Results    '''
X_opt = np.reshape(z_opt[: (N+1)*nx], (nx, N+1), order='F')  # (nx, N+1)
U_off = (N+1)*nx
U_opt = np.reshape(z_opt[U_off: U_off + N*nu], (nu, N), order='F')  # (nu, N)
T_opt = z_opt[idxT]

tX = np.linspace(0.0, T_opt, N+1)
tU = np.linspace(0.0, T_opt, N)
print("Optimal time:", T_opt)
print("Final State: ", X_opt[:, -1])

# Plot Results
fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

# positions
axs[0].plot(X_opt[0, :], X_opt[1, :])
axs[0].set_ylabel('Y [m]')
axs[0].set_xlabel('X [m]')
axs[0].grid(True); 

# velocities
axs[1].plot(tX, X_opt[2, :], label='v_x')
axs[1].plot(tX, X_opt[3, :], label='v_y')
axs[1].set_ylabel('Velocity')
axs[1].grid(True); axs[1].legend()

# controls
axs[2].step(tU, U_opt[0, :], where='post', label='u1')
axs[2].step(tU, U_opt[1, :], where='post', label='u2')
axs[2].set_ylabel('Control'); axs[2].set_xlabel('Time [s]')
axs[2].grid(True); axs[2].legend()

plt.tight_layout()
plt.show()

data = {'Xnom': X_opt, 'Unom': U_opt, 'Tnom': T_opt}
scipy.io.savemat('Guess_Trajectory.mat', data)
