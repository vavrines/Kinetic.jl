# Parallel Computing

Julia supports different categories of parallel programming natively:
- Asynchronous "tasks", or coroutines
- Multi-threading
- Distributed computing

Kinetic integrates the the latter two with the CUDA based GPU computing.



```julia
function __init__()
    np = nworkers()
    nt = Threads.nthreads()
    if nt > 1 || np > 1
        @info "Kinetic will run with $np processors and $nt threads"
    else
        @info "Kinetic will run serially"
    end

    if has_cuda()
        @info "Kinetic will run with CUDA"
        for (i, dev) in enumerate(CUDA.devices())
            @info "$i: $(CUDA.name(dev))"
        end
        @info "Scalar operation is disabled in CUDA"
        CUDA.allowscalar(false)
    end
end
```
computational resources that will be called.

