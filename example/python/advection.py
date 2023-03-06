# %%
# Julia backend
from julia.api import Julia

#%%
jl = Julia(compiled_modules=False)
from julia import KitBase as kt
import numpy as np
import matplotlib.pyplot as plt

# %%
st = kt.Setup(space="1d1f3v", collision="fsm", interpOrder=1, boundary="period", maxTime=1.0)
ps = kt.PSpace1D(0.0, 1.0, 100, 1)
vs = kt.VSpace3D(-8, 8, 48, -8, 8, 28, -8, 8, 28)

# %%
knudsen = 1e-5
muref = kt.ref_vhs_vis(knudsen, 1.0, 0.5)
fsm = kt.fsm_kernel(vs, muref, 5, 1.0)
gas = kt.Gas(Kn=knudsen, K=0.0, fsm=fsm)
phi, psi, chi = kt.kernel_mode(
    5,
    vs.u1,
    vs.v1,
    vs.w1,
    vs.du[1, 1, 1],
    vs.dv[1, 1, 1],
    vs.dw[1, 1, 1],
    vs.nu,
    vs.nv,
    vs.nw,
    1.0,
)

# %%
def fw(x, p):
    ρ = 1 + 0.1 * np.sin(2 * np.pi * x)
    u = 1.0
    λ = ρ
    return kt.prim_conserve([ρ, u, 0, 0, λ], gas.γ)

def ff(x, p):
    w = fw(x, p)
    prim = kt.conserve_prim(w, gas.γ)
    return kt.maxwellian(vs.u, vs.v, vs.w, prim)

# %%
w = np.zeros((100, 5))
f = np.zeros((100, 48, 28, 28))
df = np.zeros((100, 48, 28, 28))
Q = np.zeros((100, 48, 28, 28))
for i in range(100):
    f[i, :, :, :] = ff(ps.x[i+1], None)
    w[i, :] = kt.moments_conserve(f[i, :, :, :], vs.u, vs.v, vs.w, vs.weights)

# %%
import copy
w0 = copy.deepcopy(w)

# %%
def compute_df(df, f, ps, vs):
    for i in range(1, ps.nx - 1):
        for j in range(1, vs.nu):
            if vs.u[j, 1, 1] > 0.0:
                df[i, j, :, :] = (f[i, j, :, :] - f[i - 1, j, :, :]) / ps.dx[i]
            else:
                df[i, j, :, :] = (f[i + 1, j, :, :] - f[i, j, :, :]) / ps.dx[i]
    for j in range(1, vs.nu):
        if vs.u[j, 1, 1] > 0.0:
            df[0, j, :, :] = (f[0, j, :, :] - f[ps.nx - 1, j, :, :]) / ps.dx[0]
            df[ps.nx - 1, j, :, :] = (f[ps.nx - 1, j, :, :] - f[ps.nx - 2, j, :, :]) / ps.dx[ps.nx - 1]
        else:
            df[0, j, :, :] = (f[1, j, :, :] - f[0, j, :, :]) / ps.dx[0]
            df[ps.nx - 1, j, :, :] = (f[0, j, :, :] - f[ps.nx - 1, j, :, :]) / ps.dx[ps.nx - 1]

# %%
def compute_Q(Q, f, ps, gas, phi, psi, chi, dt):
    # f = (x, u, v, w)
    for i in range(0, ps.nx):
        Q[i, :, :, :] = dt * kt.boltzmann_fft(f[i, :, :, :], gas.fsm.Kn, gas.fsm.nm, phi, psi, chi)

def compute_Qbgk(Q, f, ps, vs, gas, muref, dt):
    for i in range(0, ps.nx):
        ww = kt.moments_conserve(f[i, :, :, :], vs.u, vs.v, vs.w, vs.weights)
        prim = kt.conserve_prim(ww, gas.γ)
        mm = kt.maxwellian(vs.u, vs.v, vs.w, prim)
        tau = kt.vhs_collision_time(prim, muref, gas.ω)
        Q[i, :, :, :] = dt * (mm - f[i, :, :, :]) / tau

# %%
def step(f, df, Q, ps, vs, dt):
    for i in range(0, ps.nx):
        f[i, :, :, :] -= vs.u * df[i, :, :, :] * dt - Q[i, :, :, :]
        w[i, :] = kt.moments_conserve(f[i, :, :, :], vs.u, vs.v, vs.w, vs.weights)

# %%
dt = 1e-3
for iter in range(15):
    compute_df(df, f, ps, vs)
    compute_Q(Q, f, ps, gas, phi, psi, chi, dt)
    #compute_Qbgk(Q, f, ps, vs, gas, muref, dt)
    step(f, df, Q, ps, vs, dt)
    print(iter + 1)

# %%
plt.plot(w[:, 1])
plt.plot(w0[:, 1], '--')
plt.legend(['w', 'w0'])
plt.show()
