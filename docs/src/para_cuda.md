# GPU computing

```julia
CUDA
```

architecture


In the following, we present an illustrative test of kinetic flux vector splitting method to evaluate upwind flux of the Boltzmann equation.
The test is conducted on a Tesla K80 GPU on [nextjournal.com](nextjournal.com).
We first load all the modules, and do a CPU-based computation.

```julia
import Pkg
Pkg.add("Revise")
Pkg.add("KitBase")
Pkg.add("CUDA")
Pkg.add("BenchmarkTools")

using Revise, CUDA, BenchmarkTools, KitBase

dt = 1e-3

primL = [1.0, 0.0, 0.5]
primR = [0.125, 0.0, 0.625]

u = collect(-5.0:0.01:5.0)
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The benchmark result on a Intel NUC8i7BEH desktop is around `5.244 μs (3 allocations: 24.00 KiB)`.
Then let's turn to GPU.
```julia
u = collect(-5.0:0.01:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The benchmark result is around `32.965 μs (187 allocations: 10.73 KiB)`.
As can be seen, due to the relative small input size, the GPU threads aren't fully occupied, and therefore CPU is more efficient in this case.

Then let's increase the input vector size, i.e. to consider more discrete particle velocity points for distribution functions.
```julia
u = collect(-5.0:0.001:5.0)
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)

u = collect(-5.0:0.001:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The results become around `50.011 μs (6 allocations: 234.80 KiB)` for CPU and `33.640 μs (187 allocations: 10.73 KiB)` for GPU.

We could further increase the computation size.
```julia
u = collect(-5.0:0.0001:5.0)
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)

u = collect(-5.0:0.0001:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@btime flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```
The results become around `507.960 μs (6 allocations: 2.29 MiB)` for CPU and `32.021 μs (187 allocations: 10.73 KiB)` for GPU.
Under this size of computation, the GPU brings about 16x efficiency increment.