# Benchmark

Here we provide a benchmark to identity the performance variation between Julia and Fortran implementations.
For brevity, we direct make use of the dynamic library [kitmod.so](https://github.com/vavrines/KitFort.jl/tree/main/src/fortran) by `ccall` function in Julia, and
compare the efficiency of computing numerical fluxes by [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).

```julia
using Kinetic, BenchmarkTools

begin
    u = collect(-5.0:0.1:5.0)
    nu = length(u)
    weights = ones(nu) .* 0.5

    fw = zeros(3)
    fh = zeros(nu)
    fb = zeros(nu)

    inK = 2
    γ = 5.0 / 3.0
    primL = [1., 0., 1.]
    wL = prim_conserve(primL, γ)
    hL = maxwellian(u, primL) |> Array;
    bL = hL .* 2 ./ (2.)
    shL = zeros(nu)
    sbL = zeros(nu)
    lenL = 0.1

    primR = [0.5, 0., 1.]
    wR = prim_conserve(primR, γ)
    hR = maxwellian(u, primR) |> Array;
    bR = hR .* 2 ./ (2.)
    shR = zeros(nu)
    sbR = zeros(nu)
    lenR = 0.1

    muref = 0.001
    omega = 0.72
    prandtl = 1.0
    dt = 1e-4
end

#--- kfvs ---#
@btime ccall(
    (:__kinetic_MOD_flux_kfvs_2f1v, "kitmod.so"),
    Nothing,
    (
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64}, 
        Ref{Int}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64},
    ),
    fw,
    fh,
    fb,
    hL,
    bL,
    hR,
    bR,
    u,
    weights,
    nu,
    dt,
    shL,
    sbL,
    shR,
    sbR,
)

@btime flux_kfvs!(fw, fh, fb, hL, bL, hR, bR, u, weights, dt, shL, sbL, shR, sbR)

#--- ugks ---#
@btime ccall(
    (:__kinetic_MOD_flux_ugks_2f1v, "kitmod.so"),
    Nothing,
    (
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Int}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64}, 
        Ref{Float64},
    ),
    fw,
    fh,
    fb,
    wL,
    hL,
    bL,
    wR,
    hR,
    bR,
    u,
    weights,
    nu,
    inK,
    γ,
    muref,
    omega,
    prandtl,
    dt,
    lenL,
    lenR,
    shL,
    sbL,
    shR,
    sbR,
)

@btime flux_ugks!(fw, fh, fb, wL, hL, bL, wR, hR, bR, u, weights, inK, γ, muref, omega, prandtl, dt, lenL, lenR, shL, sbL, shR, sbR)

```

The results on a intel NUC8i7BEH with i7-8559U with 101 velocity points is as follows

Kinetic.jl
- KFVS flux ~ 6.747 μs (13 allocations: 11.38 KiB)
- UGKS flux ~ 13.344 μs (123 allocations: 20.94 KiB)

KitFort.jl
- KFVS flux ~ 5.421 μs (37 allocations: 800 bytes)
- UGKS flux ~ 11.413 μs (55 allocations: 1.09 KiB)

As presented, there is an improvement on efficiency by around 15%.