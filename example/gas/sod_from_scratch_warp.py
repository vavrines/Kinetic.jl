# %%
"""
Sod shock tube for the 1D reduced BGK equation, rewritten with NVIDIA Warp.
This follows the structure of sod_from_scratch.jl and the NumPy translation sod_from_scratch.py.

The main difference is data layout. A list of small Python objects is not a good
GPU data structure, so each field is stored as one contiguous Warp array:

  ctr.w     : shape nx*3      flattened as w[i*3 + k]
  ctr.prim  : shape nx*3      flattened as prim[i*3 + k]
  ctr.h,b   : shape nx*nu     flattened as h[i*nu + j]
  face.fw   : shape (nx+1)*3  flattened as fw[i*3 + k]
  face.fh,fb: shape (nx+1)*nu flattened as fh[i*nu + j]

Run as a script:
  python sod_from_scratch_warp.py

Optional script arguments:
  python sod_from_scratch_warp.py --device cuda
  python sod_from_scratch_warp.py --device cpu

Run in a Jupyter notebook or VSCode cell:
  x, sol, ks, ctr, face = run_sod(device="cuda", nx=200, nu=72, show=True)
  x, sol, ks, ctr, face = run_sod(device="cpu", nx=200, nu=72, show=True)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import warp as wp


# --- Data Classes --- #


@dataclass
class Setup:
    # Solver configuration
    cfl: float
    maxTime: float


@dataclass
class PSpace1D:
    # Physical mesh
    x0: float
    x1: float
    nx: int
    x: wp.array
    dx: wp.array


def PSpace1D_factory(x0: float, x1: float, nx: int, device: str) -> PSpace1D:
    # Build a uniform cell-centered grid in physical space on host,
    # then move arrays to the Warp device.
    delta = (x1 - x0) / nx
    x_np = np.empty(nx, dtype=np.float64)
    dx_np = np.empty(nx, dtype=np.float64)
    for i in range(nx):
        x_np[i] = x0 + (i + 0.5) * delta
        dx_np[i] = delta
    return PSpace1D(
        x0,
        x1,
        nx,
        wp.array(x_np, dtype=wp.float64, device=device),
        wp.array(dx_np, dtype=wp.float64, device=device),
    )


@dataclass
class VSpace1D:
    # Velocity mesh
    u0: float
    u1: float
    nu: int
    u: wp.array
    du: wp.array
    weights: wp.array


def VSpace1D_factory(u0: float, u1: float, nu: int, device: str) -> VSpace1D:
    # Build a uniform quadrature grid in velocity space.
    delta = (u1 - u0) / nu
    u_np = np.empty(nu, dtype=np.float64)
    du_np = np.empty(nu, dtype=np.float64)
    weights_np = np.empty(nu, dtype=np.float64)
    for i in range(nu):
        u_np[i] = u0 + (i + 0.5) * delta
        du_np[i] = delta
        weights_np[i] = delta
    return VSpace1D(
        u0,
        u1,
        nu,
        wp.array(u_np, dtype=wp.float64, device=device),
        wp.array(du_np, dtype=wp.float64, device=device),
        wp.array(weights_np, dtype=wp.float64, device=device),
    )


@dataclass
class Gas:
    # Gas property
    Kn: float
    K: float = 2.0
    gamma: float = 5.0 / 3.0
    omega: float = 0.81
    alpha_r: float = 1.0
    omega_r: float = 0.5
    mu_r: Optional[float] = None

    def __post_init__(self) -> None:
        if self.mu_r is None:
            self.mu_r = ref_vhs_vis(self.Kn, self.alpha_r, self.omega_r)


@dataclass
class ControlVolume2F:
    # Macroscopic state and mesoscopic distributions in all cells.
    w: wp.array  # nx * 3
    prim: wp.array  # nx * 3
    h: wp.array  # nx * nu
    b: wp.array  # nx * nu


@dataclass
class Interface2F:
    # Fluxes across all cell interfaces.
    fw: wp.array  # (nx + 1) * 3
    fh: wp.array  # (nx + 1) * nu
    fb: wp.array  # (nx + 1) * nu


@dataclass
class SolverSet:
    # Top-level solver state container
    set: Setup
    ps: PSpace1D
    vs: VSpace1D
    gas: Gas
    device: str


# --- Host-side scalar physics --- #


def ref_vhs_vis(Kn: float, alpha: float, omega: float) -> float:
    # Reference viscosity model for variable hard sphere gas.
    return (
        5.0
        * (alpha + 1.0)
        * (alpha + 2.0)
        * math.sqrt(math.pi)
        / (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega))
        * Kn
    )


# --- Warp scalar functions --- #


@wp.func
def _prim_conserve_0(
    rho: wp.float64, U: wp.float64, lam: wp.float64, gamma: wp.float64
) -> wp.float64:
    return rho


@wp.func
def _prim_conserve_1(
    rho: wp.float64, U: wp.float64, lam: wp.float64, gamma: wp.float64
) -> wp.float64:
    return rho * U


@wp.func
def _prim_conserve_2(
    rho: wp.float64, U: wp.float64, lam: wp.float64, gamma: wp.float64
) -> wp.float64:
    return (
        wp.float64(0.5) * rho / lam / (gamma - wp.float64(1.0))
        + wp.float64(0.5) * rho * U * U
    )


@wp.func
def _conserve_prim_0(
    w0: wp.float64, w1: wp.float64, w2: wp.float64, gamma: wp.float64
) -> wp.float64:
    return w0


@wp.func
def _conserve_prim_1(
    w0: wp.float64, w1: wp.float64, w2: wp.float64, gamma: wp.float64
) -> wp.float64:
    return w1 / w0


@wp.func
def _conserve_prim_2(
    w0: wp.float64, w1: wp.float64, w2: wp.float64, gamma: wp.float64
) -> wp.float64:
    return (
        wp.float64(0.5)
        * w0
        / (gamma - wp.float64(1.0))
        / (w2 - wp.float64(0.5) * w1 * w1 / w0)
    )


@wp.func
def _maxwellian_h(
    u: wp.float64, rho: wp.float64, U: wp.float64, lam: wp.float64
) -> wp.float64:
    return rho * wp.sqrt(lam / wp.float64(wp.pi)) * wp.exp(-lam * (u - U) * (u - U))


@wp.func
def _maxwellian_b(h: wp.float64, K: wp.float64, lam: wp.float64) -> wp.float64:
    return h * K / (wp.float64(2.0) * lam)


@wp.func
def _vhs_collision_time(
    rho: wp.float64, lam: wp.float64, mu_r: wp.float64, omega: wp.float64
) -> wp.float64:
    return mu_r * wp.float64(2.0) * wp.pow(lam, wp.float64(1.0) - omega) / rho


@wp.func
def _sound_speed(lam: wp.float64, gamma: wp.float64) -> wp.float64:
    return wp.sqrt(wp.float64(0.5) * gamma / lam)


@wp.func
def _ic_sod_rho(x: wp.float64) -> wp.float64:
    if x < wp.float64(0.5):
        return wp.float64(1.0)
    return wp.float64(0.125)


@wp.func
def _ic_sod_u(x: wp.float64) -> wp.float64:
    return wp.float64(0.0)


@wp.func
def _ic_sod_lambda(x: wp.float64) -> wp.float64:
    if x < wp.float64(0.5):
        return wp.float64(0.5)
    return wp.float64(0.625)


# --- Kernels --- #


@wp.kernel
def _init_fvm_kernel(
    x: wp.array(dtype=wp.float64),
    u: wp.array(dtype=wp.float64),
    w: wp.array(dtype=wp.float64),
    prim: wp.array(dtype=wp.float64),
    h: wp.array(dtype=wp.float64),
    b: wp.array(dtype=wp.float64),
    nx: int,
    nu: int,
    K: wp.float64,
    gamma: wp.float64,
):
    i, j = wp.tid()

    xi = x[i]
    rho = _ic_sod_rho(xi)
    U = _ic_sod_u(xi)
    lam = _ic_sod_lambda(xi)

    if j == 0:
        base3 = i * 3
        prim[base3 + 0] = rho
        prim[base3 + 1] = U
        prim[base3 + 2] = lam
        w[base3 + 0] = _prim_conserve_0(rho, U, lam, gamma)
        w[base3 + 1] = _prim_conserve_1(rho, U, lam, gamma)
        w[base3 + 2] = _prim_conserve_2(rho, U, lam, gamma)

    uj = u[j]
    mh = _maxwellian_h(uj, rho, U, lam)
    mb = _maxwellian_b(mh, K, lam)
    base = i * nu + j
    h[base] = mh
    b[base] = mb


@wp.kernel
def _local_speed_kernel(
    prim: wp.array(dtype=wp.float64),
    dx: wp.array(dtype=wp.float64),
    local_speed: wp.array(dtype=wp.float64),
    nx: int,
    u1: wp.float64,
    gamma: wp.float64,
):
    i = wp.tid()
    base3 = i * 3
    U = prim[base3 + 1]
    lam = prim[base3 + 2]
    sos = _sound_speed(lam, gamma)
    vmax = wp.max(u1, wp.abs(U)) + sos
    local_speed[i] = vmax / dx[i]


@wp.kernel
def _evolve_kernel(
    fw: wp.array(dtype=wp.float64),
    fh: wp.array(dtype=wp.float64),
    fb: wp.array(dtype=wp.float64),
    h: wp.array(dtype=wp.float64),
    b: wp.array(dtype=wp.float64),
    u: wp.array(dtype=wp.float64),
    weights: wp.array(dtype=wp.float64),
    nx: int,
    nu: int,
    dt: wp.float64,
):
    # Each thread owns one interior interface i = 1, ..., nx-1.
    # This deliberately mirrors evolve! in the original script.
    tid = wp.tid()
    iface = tid + 1

    sum0 = wp.float64(0.0)
    sum1 = wp.float64(0.0)
    sum2h = wp.float64(0.0)
    sum2b = wp.float64(0.0)

    left_cell = iface - 1
    right_cell = iface

    for j in range(nu):
        uj = u[j]
        wt = weights[j]
        hL = h[left_cell * nu + j]
        bL = b[left_cell * nu + j]
        hR = h[right_cell * nu + j]
        bR = b[right_cell * nu + j]

        if uj >= wp.float64(0.0):
            hup = hL
            bup = bL
        else:
            hup = hR
            bup = bR

        fh[iface * nu + j] = dt * uj * hup
        fb[iface * nu + j] = dt * uj * bup

        sum0 += wt * uj * hup
        sum1 += wt * uj * uj * hup
        sum2h += wt * uj * uj * uj * hup
        sum2b += wt * uj * bup

    base3 = iface * 3
    fw[base3 + 0] = dt * sum0
    fw[base3 + 1] = dt * sum1
    fw[base3 + 2] = dt * wp.float64(0.5) * (sum2h + sum2b)


@wp.kernel
def _update_kernel(
    w: wp.array(dtype=wp.float64),
    prim: wp.array(dtype=wp.float64),
    h: wp.array(dtype=wp.float64),
    b: wp.array(dtype=wp.float64),
    fw: wp.array(dtype=wp.float64),
    fh: wp.array(dtype=wp.float64),
    fb: wp.array(dtype=wp.float64),
    u: wp.array(dtype=wp.float64),
    dx: wp.array(dtype=wp.float64),
    nx: int,
    nu: int,
    K: wp.float64,
    gamma: wp.float64,
    mu_r: wp.float64,
    omega: wp.float64,
    dt: wp.float64,
):
    # Each thread owns one interior cell i = 1, ..., nx-2.
    # Boundary cells are kept fixed, as in the original update! loop.
    tid = wp.tid()
    i = tid + 1

    base3 = i * 3
    left3 = i * 3
    right3 = (i + 1) * 3
    dxi = dx[i]

    w0 = w[base3 + 0] + (fw[left3 + 0] - fw[right3 + 0]) / dxi
    w1 = w[base3 + 1] + (fw[left3 + 1] - fw[right3 + 1]) / dxi
    w2 = w[base3 + 2] + (fw[left3 + 2] - fw[right3 + 2]) / dxi

    w[base3 + 0] = w0
    w[base3 + 1] = w1
    w[base3 + 2] = w2

    rho = _conserve_prim_0(w0, w1, w2, gamma)
    U = _conserve_prim_1(w0, w1, w2, gamma)
    lam = _conserve_prim_2(w0, w1, w2, gamma)

    prim[base3 + 0] = rho
    prim[base3 + 1] = U
    prim[base3 + 2] = lam

    tau = _vhs_collision_time(rho, lam, mu_r, omega)
    factor = dt / tau
    denom = wp.float64(1.0) + factor

    for j in range(nu):
        uj = u[j]
        MH = _maxwellian_h(uj, rho, U, lam)
        MB = _maxwellian_b(MH, K, lam)

        cidx = i * nu + j
        lidx = i * nu + j
        ridx = (i + 1) * nu + j

        h[cidx] = (h[cidx] + (fh[lidx] - fh[ridx]) / dxi + factor * MH) / denom
        b[cidx] = (b[cidx] + (fb[lidx] - fb[ridx]) / dxi + factor * MB) / denom


# --- Solver functions with the same high-level names as the source scripts --- #


def init_fvm(ks: SolverSet) -> tuple[ControlVolume2F, Interface2F]:
    nx = ks.ps.nx
    nu = ks.vs.nu
    device = ks.device

    ctr = ControlVolume2F(
        w=wp.zeros(nx * 3, dtype=wp.float64, device=device),
        prim=wp.zeros(nx * 3, dtype=wp.float64, device=device),
        h=wp.zeros(nx * nu, dtype=wp.float64, device=device),
        b=wp.zeros(nx * nu, dtype=wp.float64, device=device),
    )
    face = Interface2F(
        fw=wp.zeros((nx + 1) * 3, dtype=wp.float64, device=device),
        fh=wp.zeros((nx + 1) * nu, dtype=wp.float64, device=device),
        fb=wp.zeros((nx + 1) * nu, dtype=wp.float64, device=device),
    )

    wp.launch(
        _init_fvm_kernel,
        dim=(nx, nu),
        inputs=[
            ks.ps.x,
            ks.vs.u,
            ctr.w,
            ctr.prim,
            ctr.h,
            ctr.b,
            nx,
            nu,
            ks.gas.K,
            ks.gas.gamma,
        ],
        device=device,
    )
    return ctr, face


def timestep(ks: SolverSet, ctr: ControlVolume2F) -> float:
    nx = ks.ps.nx
    local_speed = wp.zeros(nx, dtype=wp.float64, device=ks.device)
    wp.launch(
        _local_speed_kernel,
        dim=nx,
        inputs=[ctr.prim, ks.ps.dx, local_speed, nx, ks.vs.u1, ks.gas.gamma],
        device=ks.device,
    )
    # Keeping the global max on host keeps this example simple and transparent.
    # For production, replace this with a Warp reduction.
    tmax = float(np.max(local_speed.numpy()))
    return ks.set.cfl / tmax


def evolve(ks: SolverSet, ctr: ControlVolume2F, face: Interface2F, dt: float) -> None:
    wp.launch(
        _evolve_kernel,
        dim=ks.ps.nx - 1,
        inputs=[
            face.fw,
            face.fh,
            face.fb,
            ctr.h,
            ctr.b,
            ks.vs.u,
            ks.vs.weights,
            ks.ps.nx,
            ks.vs.nu,
            dt,
        ],
        device=ks.device,
    )


def update(ks: SolverSet, ctr: ControlVolume2F, face: Interface2F, dt: float) -> None:
    wp.launch(
        _update_kernel,
        dim=ks.ps.nx - 2,
        inputs=[
            ctr.w,
            ctr.prim,
            ctr.h,
            ctr.b,
            face.fw,
            face.fh,
            face.fb,
            ks.vs.u,
            ks.ps.dx,
            ks.ps.nx,
            ks.vs.nu,
            ks.gas.K,
            ks.gas.gamma,
            ks.gas.mu_r,
            ks.gas.omega,
            dt,
        ],
        device=ks.device,
    )


def extract_sol(ps: PSpace1D, ctr: ControlVolume2F) -> np.ndarray:
    # Gather primitive variables from the device into an nx-by-3 host array.
    return ctr.prim.numpy().reshape(ps.nx, 3).copy()


def run_sod(
    device: str = "cuda",
    nx: int = 200,
    nu: int = 72,
    cfl: float = 0.5,
    max_time: float = 0.2,
    x0: float = 0.0,
    x1: float = 1.0,
    u0: float = -5.0,
    u1: float = 5.0,
    kn: float = 1.0e-4,
    show: bool = True,
    return_state: bool = True,
):
    """Run the Sod example without using argparse.

    This is the preferred entry point for Jupyter notebooks and VSCode cells.

    Example
    -------
    x, sol, ks, ctr, face = run_sod(device="cuda", nx=200, nu=72, show=True)

    Returns
    -------
    If return_state=True:
        x, sol, ks, ctr, face
    otherwise:
        x, sol
    where sol[:, 0], sol[:, 1], and 0.5*sol[:,0]/sol[:,2] are density,
    velocity, and pressure, respectively.
    """
    # Calling wp.init() repeatedly is safe in interactive sessions.
    wp.init()

    # Set up physical, velocity, and gas parameters.
    set_ = Setup(cfl, max_time)
    ps = PSpace1D_factory(x0, x1, nx, device)
    vs = VSpace1D_factory(u0, u1, nu, device)
    gas = Gas(Kn=kn)
    ks = SolverSet(set_, ps, vs, gas, device)

    # Initialize the finite-volume state and face flux arrays.
    ctr, face = init_fvm(ks)

    # Compute a stable time step and number of iterations.
    dt = timestep(ks, ctr)
    nt = int(ks.set.maxTime // dt)

    # Main time-stepping loop.
    for _ in range(nt):
        evolve(ks, ctr, face, dt)
        update(ks, ctr, face, dt)

    wp.synchronize_device(device)

    # Extract primitive solution.
    sol = extract_sol(ps, ctr)
    x = ps.x.numpy()

    if show:
        plt.figure()
        plt.plot(x, sol[:, 0], label="density", linewidth=1.5)
        plt.plot(x, sol[:, 1], label="velocity", linewidth=1.5)
        plt.plot(x, 0.5 * sol[:, 0] / sol[:, 2], label="pressure", linewidth=1.5)
        plt.xlabel("x")
        plt.legend()
        plt.tight_layout()
        plt.show()

    if return_state:
        return x, sol, ks, ctr, face
    return x, sol


def main(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--device", default="cuda", help="Warp device, e.g. cuda, cuda:0, or cpu"
    )
    parser.add_argument("--nx", type=int, default=200)
    parser.add_argument("--nu", type=int, default=72)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)

    run_sod(
        device=args.device,
        nx=args.nx,
        nu=args.nu,
        show=not args.no_show,
        return_state=False,
    )


# %%
if __name__ == "__main__":
    # main()
    x, sol, ks, ctr, face = run_sod(device="cuda", nx=200, nu=72, show=True)
