using KitBase, CUDA

u = collect(-5.0:0.0001:5.0)
primL = [1.0, 0.0, 0.5]
primR = [0.125, 0.0, 0.625]
dt = 1.0

fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@time flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```0.001932 seconds (6 allocations: 2.289 MiB)```

u = collect(-5.0:0.0001:5.0) |> CuArray
fL = maxwellian(u, primL)
fR = maxwellian(u, primR)
ff = similar(fL)
sfL = zero(fL)
sfR = zero(fR)

@time flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)
```0.000219 seconds (203 allocations: 15.125 KiB)```
