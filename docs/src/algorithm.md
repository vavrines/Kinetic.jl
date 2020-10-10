# Solution Algorithm

The finite volume method (FVM) is employed in this package. The general solution algorithm can be conclude as follows, where both explicit and implicit methods are implemented.

![](./assets/solver_process.png)

(1) pre-process: `initialize(configfilename::AbstractString)`
* new run: .txt / .cfg / .toml
* restart: .jld2 (HDF5) -> solSet, ctr, face, t

(2) calculate timestep based on the CFL condition

(3) reconstruct field solutions
* van-Leer
* minmod
* superbee
* WENO

(4) evolve numerical fluxes
* macroscopic: Godunov, Lax, Roe, HLL, wave-propagation
* mesoscopic: upwind, central-upwind, gas-kinetic scheme

(5) update cell-averaged variables
* explicit
* implicit
* implicit-explicit (IMEX)

(6) post-process