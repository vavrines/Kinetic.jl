# Multiple threading

The multi-threading computation is built upon Julia's `@threads` macro.
```julia
Base.Threads.@threads for ... end
```
It provides an OpenMP type parallelization.
The iteration space is splitted among multiple tasks and those tasks are parallelized on threads according to a scheduling policy.
A barrier is placed at the end of the loop which waits for all tasks to finish execution.

In Kinetic, `@threads` is set in front of the loops for reconstruction, evolution and update.
For example, the evaluation of fluxes is conducted as follows.
```julia
@inbounds Threads.@threads for i = idx0:idx1
    flux_gks!(
        face[i].fw,
        ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw,
        ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw,
        KS.gas.γ,
        KS.gas.K,
        KS.gas.μᵣ,
        KS.gas.ω,
        dt,
        0.5 * ctr[i-1].dx,
        0.5 * ctr[i].dx,
        ctr[i-1].sw,
        ctr[i].sw,
    )
end
```
It automatically makes use of multiple threading if Julia is initialized with
```bash
julia -t n
```

Besides of `@threads`, finer dispatch can be made with `@spawn` macro.
```julia
Base.Threads.@spawn
```
It creates and runs a task on any available thread.
To wait for the task to finish, call `wait` on the result of this macro, or call fetch to wait and then obtain its return value.
Values can be interpolated into `@spawn` via `$`, which copies the value directly into the constructed underlying closure.
It allows user to insert the value of a variable, isolating the aysnchronous code from changes to the variable's value in the current task.
This can be conduct with the low-level reconstruction, flux and step functions.