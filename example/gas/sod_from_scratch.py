# %%
"""
Translated Sod shock tube example from Julia to Python using NumPy.
This script retains the original structure from `sod_from_scratch.jl`.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Optional

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
    x: np.ndarray
    dx: np.ndarray


def PSpace1D_factory(x0: float, x1: float, nx: int) -> PSpace1D:
    # Build a uniform cell-centered grid in physical space
    delta = (x1 - x0) / nx
    x = np.empty(nx, dtype=float)
    dx = np.empty(nx, dtype=float)

    for i in range(nx):
        x[i] = x0 + (i + 0.5) * delta
        dx[i] = delta

    return PSpace1D(x0, x1, nx, x, dx)


@dataclass
class VSpace1D:
    # Velocity mesh
    u0: float
    u1: float
    nu: int
    u: np.ndarray
    du: np.ndarray
    weights: np.ndarray


def VSpace1D_factory(u0: float, u1: float, nu: int) -> VSpace1D:
    # Build a uniform quadrature grid in velocity space
    delta = (u1 - u0) / nu
    u = np.empty(nu, dtype=float)
    du = np.empty(nu, dtype=float)
    weights = np.empty(nu, dtype=float)

    for i in range(nu):
        u[i] = u0 + (i + 0.5) * delta
        du[i] = delta
        weights[i] = delta

    return VSpace1D(u0, u1, nu, u, du, weights)


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

    def __post_init__(self):
        # If the reference viscosity was not provided, compute it from model constants
        if self.mu_r is None:
            self.mu_r = ref_vhs_vis(self.Kn, self.alpha_r, self.omega_r)


@dataclass
class ControlVolume2F:
    # Macroscopic state and mesoscopic distributions in each cell
    w: np.ndarray
    prim: np.ndarray
    h: np.ndarray
    b: np.ndarray


@dataclass
class Interface2F:
    # Fluxes across each cell interface
    fw: np.ndarray
    fh: np.ndarray
    fb: np.ndarray


@dataclass
class SolverSet:
    # Top-level solver state container
    set: Setup
    ps: PSpace1D
    vs: VSpace1D
    gas: Gas


# --- Physics --- #


def conserve_prim(W: np.ndarray, gamma: float) -> np.ndarray:
    # Convert conservative variables (rho, rho*u, energy) to primitive form (rho, u, lambda)
    prim = np.zeros_like(W)
    prim[0] = W[0]
    prim[1] = W[1] / W[0]
    prim[2] = 0.5 * W[0] / (gamma - 1.0) / (W[2] - 0.5 * W[1] ** 2 / W[0])
    return prim


def prim_conserve(prim: np.ndarray, gamma: float) -> np.ndarray:
    # Convert primitive variables back to conserved variables
    W = np.zeros_like(prim)
    W[0] = prim[0]
    W[1] = prim[0] * prim[1]
    W[2] = 0.5 * prim[0] / prim[2] / (gamma - 1.0) + 0.5 * prim[0] * prim[1] ** 2
    return W


def maxwellian(u: np.ndarray, prim: np.ndarray, K: float):
    # Reduced Maxwellian distribution and related moment term
    rho, U, lam = prim
    h = rho * np.sqrt(lam / np.pi) * np.exp(-lam * (u - U) ** 2)
    b = h * K / (2.0 * lam)
    return h, b


def ref_vhs_vis(Kn: float, alpha: float, omega: float) -> float:
    # Reference viscosity model for variable hard sphere gas
    return (
        5.0
        * (alpha + 1.0)
        * (alpha + 2.0)
        * np.sqrt(np.pi)
        / (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega))
        * Kn
    )


def vhs_collision_time(prim: np.ndarray, mu_r: float, omega: float) -> float:
    # Collision time in the VHS kinetic model
    rho = prim[0]
    lam = prim[-1]
    return mu_r * 2.0 * lam ** (1.0 - omega) / rho


def sound_speed(prim: np.ndarray, gamma: float) -> float:
    # Sound speed from primitive variables
    return np.sqrt(0.5 * gamma / prim[-1])


# --- Flux --- #


def heaviside(arr: np.ndarray) -> np.ndarray:
    # Elementwise step function used to split flux contributions by velocity sign
    return (arr >= 0.0).astype(float)


def flux_kfvs(
    fw: np.ndarray,
    fh: np.ndarray,
    fb: np.ndarray,
    hL: np.ndarray,
    bL: np.ndarray,
    hR: np.ndarray,
    bR: np.ndarray,
    u: np.ndarray,
    weights: np.ndarray,
    dt: float,
) -> None:
    # Kinetic flux vector splitting: choose upstream distribution for each velocity
    delta = heaviside(u)
    h = hL * delta + hR * (1.0 - delta)
    b = bL * delta + bR * (1.0 - delta)

    # Update mesoscopic fluxes in velocity space
    fh[:] = dt * u * h
    fb[:] = dt * u * b

    # Compute macroscopic flux moments
    fw[0] = dt * np.sum(weights * u * h)
    fw[1] = dt * np.sum(weights * u ** 2 * h)
    fw[2] = dt * 0.5 * (np.sum(weights * u ** 3 * h) + np.sum(weights * u * b))


def evolve(
    ks: SolverSet, ctr: List[ControlVolume2F], face: List[Interface2F], dt: float
) -> None:
    # Compute fluxes at each cell interface from neighboring cell states
    for i in range(1, ks.ps.nx):
        flux_kfvs(
            face[i].fw,
            face[i].fh,
            face[i].fb,
            ctr[i - 1].h,
            ctr[i - 1].b,
            ctr[i].h,
            ctr[i].b,
            ks.vs.u,
            ks.vs.weights,
            dt,
        )


# --- Update --- #


def timestep(ks: SolverSet, ctr: List[ControlVolume2F]) -> float:
    # Compute stable time step using CFL condition and local wave speeds
    tmax = 0.0
    for i in range(ks.ps.nx):
        prim = ctr[i].prim
        sos = sound_speed(prim, ks.gas.gamma)
        vmax = max(ks.vs.u1, abs(prim[1])) + sos
        tmax = max(tmax, vmax / ks.ps.dx[i])
    return ks.set.cfl / tmax


def step(
    w: np.ndarray,
    prim: np.ndarray,
    h: np.ndarray,
    b: np.ndarray,
    fwL: np.ndarray,
    fhL: np.ndarray,
    fbL: np.ndarray,
    fwR: np.ndarray,
    fhR: np.ndarray,
    fbR: np.ndarray,
    u: np.ndarray,
    K: float,
    gamma: float,
    mu_r: float,
    omega: float,
    dx: float,
    dt: float,
) -> None:
    # Update conservative state using interface fluxes,
    # then relax the kinetic distributions toward the local Maxwellian.
    w += (fwL - fwR) / dx
    prim[:] = conserve_prim(w, gamma)

    MH, MB = maxwellian(u, prim, K)
    tau = vhs_collision_time(prim, mu_r, omega)

    factor = dt / tau
    denom = 1.0 + factor
    for i in range(len(u)):
        h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + factor * MH[i]) / denom
        b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + factor * MB[i]) / denom


def update(
    ks: SolverSet, ctr: List[ControlVolume2F], face: List[Interface2F], dt: float
) -> None:
    # Advance each cell's macroscopic and mesoscopic state using precomputed interface fluxes
    for i in range(1, ks.ps.nx - 1):
        step(
            ctr[i].w,
            ctr[i].prim,
            ctr[i].h,
            ctr[i].b,
            face[i].fw,
            face[i].fh,
            face[i].fb,
            face[i + 1].fw,
            face[i + 1].fh,
            face[i + 1].fb,
            ks.vs.u,
            ks.gas.K,
            ks.gas.gamma,
            ks.gas.mu_r,
            ks.gas.omega,
            ks.ps.dx[i],
            dt,
        )


# --- IO --- #


def ic_sod(x: float) -> np.ndarray:
    # Sod shock tube initial condition: left and right states separated at x=0.5
    if x < 0.5:
        return np.array([1.0, 0.0, 0.5], dtype=float)
    return np.array([0.125, 0.0, 0.625], dtype=float)


def init_fvm(ks: SolverSet):
    # Initialize cell-averaged variables and interface flux storage
    gamma = ks.gas.gamma
    K = ks.gas.K
    ps = ks.ps
    vs = ks.vs

    ctr: List[ControlVolume2F] = []
    face: List[Interface2F] = []

    for i in range(ps.nx):
        prim = ic_sod(ps.x[i])
        w = prim_conserve(prim, gamma)
        h, b = maxwellian(vs.u, prim, K)
        ctr.append(ControlVolume2F(w, prim, h, b))

    for _ in range(ps.nx + 1):
        fw = np.zeros(3, dtype=float)
        fh = np.zeros(vs.nu, dtype=float)
        fb = np.zeros(vs.nu, dtype=float)
        face.append(Interface2F(fw, fh, fb))

    return ctr, face


def extract_sol(ps: PSpace1D, ctr: List[ControlVolume2F]) -> np.ndarray:
    # Gather primitive variables from each control volume into a single array
    sol = np.zeros((ps.nx, ctr[0].prim.shape[0]), dtype=float)
    for i in range(ps.nx):
        sol[i, :] = ctr[i].prim
    return sol


def main() -> None:
    # Set up physical, velocity, and gas parameters
    set_ = Setup(0.5, 0.2)
    ps = PSpace1D_factory(0.0, 1.0, 200)
    vs = VSpace1D_factory(-5.0, 5.0, 72)
    gas = Gas(Kn=1e-4)
    ks = SolverSet(set_, ps, vs, gas)

    # Initialize the finite-volume state and face flux arrays
    ctr, face = init_fvm(ks)

    # Compute a stable time step and number of iterations
    dt = timestep(ks, ctr)
    nt = int(ks.set.maxTime // dt)

    # Main time-stepping loop
    for _ in range(nt):
        evolve(ks, ctr, face, dt)
        update(ks, ctr, face, dt)

    # Extract primitive solution and plot density/velocity/pressure
    sol = extract_sol(ps, ctr)

    plt.plot(ks.ps.x, sol[:, 0], label="density", linewidth=1.5)
    plt.plot(ks.ps.x, sol[:, 1], label="velocity", linewidth=1.5)
    plt.plot(ks.ps.x, 0.5 * sol[:, 0] / sol[:, 2], label="pressure", linewidth=1.5)
    plt.xlabel("x")
    plt.legend()
    plt.tight_layout()
    plt.show()


# %%
if __name__ == "__main__":
    main()
