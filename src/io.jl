# ============================================================
# Input Methods
# ============================================================

export read_dict, write_jld, plot_line


"""
Read text into dictionary

* @param[in]  filename   :   configuration text file
* @param[in]  allowed    :   keywords
* @return     vars       :   dictionary with values of variables

"""
function read_dict(filename::String, allowed)

    #println("Reading config from $filename")
    f = open(filename)
    vars = Dict{String,Any}()

    for line in eachline(f)
        # skip comments
        if length(line) == 0 || line[1] == '#'
            #println("skip comment line")
            continue
        end

        #print("\t");
        #println(line)
        var, val = split(line, "=")
        stripped = strip(var)
        if stripped in allowed
            #println(line)

            #vars[stripped] = parse(Float64, val)
            #vars[stripped] = strip(val)
            tmp = tryparse(Float64, val)
            if isa(tmp, Nothing)
                vars[stripped] = strip(val)
            else
                vars[stripped] = isinteger(tmp) ? Int(tmp) : tmp
            end
        end
    end

    println("")
    return vars

end


function read_dict(filename::String)

    f = open(filename)
    vars = Dict{String,Any}()

    for line in eachline(f)
        # skip comments
        if length(line) == 0 || line[1] == '#'
            continue
        end

        var, val = split(line, "=")
        stripped = strip(var)
        println(line)

        tmp = tryparse(Float64, val)
        if isa(tmp, Nothing)
            vars[stripped] = strip(val)
        else
            vars[stripped] = isinteger(tmp) ? Int(tmp) : tmp
        end
    end

    println("")
    return vars

end


"""
Write data into JLD2

`write_jld(KS, ctr, t)`

"""
function write_jld(
    KS::AbstractSolverSet,
    ctr::AbstractArray{<:AbstractControlVolume,1},
    t = 0::Real,
)

    strIter = string(t)
    fileOut = KS.outputFolder * "data/t=" * strIter * ".jld2"

    @save fileOut KS ctr t

end


"""
Plot 1D profile

`plot_line(KS, ctr; backend)`

"""
function plot_line(
    KS::AbstractSolverSet,
    ctr::AbstractArray{<:AbstractControlVolume1D,1};
    backend = :plots::Symbol,
)

    pltx = KS.pSpace.x[1:KS.pSpace.nx]
    plty = zeros(KS.pSpace.nx, 6)

    for i in eachindex(pltx)
        for j = 1:2
            plty[i, j] = ctr[i].prim[j]
        end

        plty[i, 3] = 1.0 / ctr[i].prim[end]
    end

    if backend == :plots
        p1 = plot(pltx, plty[:, 1], label = "Density", lw = 2, xlabel = "x")
        p1 = plot!(pltx, plty[:, 2], label = "Velocity", lw = 2)
        p1 = plot!(pltx, plty[:, 3], label = "Temperature", lw = 2)
        display(p1)
    elseif backend == :gr
        xlabel("x")
        ylabel("Density")
        legend("n")
        p1 = plot(pltx, plty[:, 1])
        display(p1)
        xlabel("x")
        ylabel("Velocity")
        legend("U")
        p2 = plot(pltx, plty[:, 2])
        display(p2)
        xlabel("x")
        ylabel("Temperature")
        legend("T")
        p3 = plot(pltx, plty[:, 3])
        display(p3)
    else
        throw("undefined plotting backend")
    end

end
