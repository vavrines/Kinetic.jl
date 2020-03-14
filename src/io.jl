# ============================================================
# Methods of Input
# ============================================================


export read_dict,
       write_jld,
       plot_line


"""
# Read text into dictionary
#
# >@param[in]  filename     :  configuration text file
# >@param[in]  allowed      :  keywords
# >@return  vars            :  dictionary with values of variables
"""
function read_dict(filename::String, allowed)

    #println("Reading config from $filename")
    f = open(filename)
    vars = Dict{String, Any}()

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
            println(line)

            #vars[stripped] = parse(Float64, val)
            vars[stripped] = strip(val)
        end
    end

    println("")
    return vars

end


# ------------------------------------------------------------
# Write output data with JLD2
# ------------------------------------------------------------
function write_jld(KS::AbstractSolverSet, ctr::AbstractArray{<:AbstractControlVolume,1}, t::Real)

    strIter = string(t)
    fileOut = KS.outputFolder * "data/t=" * strIter * ".jld2"

    @save fileOut KS ctr t

end


# ------------------------------------------------------------
# Plot line
# ------------------------------------------------------------
function plot_line(KS::AbstractSolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1})

    pltx = deepcopy(KS.pMesh.x)
    plty = zeros(KS.pMesh.nx, 6)

    for i in eachindex(pltx)
        for j=1:2
            plty[i,j] = ctr[i].prim[j]
        end

        plty[i,3] = 1. / ctr[i].prim[3]
    end

    xlabel("X"); ylabel("Density")
    #legend("N")
    p1 = plot(pltx, plty[:,1])
    display(p1)

    xlabel("X"); ylabel("Velocity")
    #legend("U")
    p2 = plot(pltx, plty[:,2])
    display(p2)

    xlabel("X"); ylabel("Temperature")
    #legend("T")
    p3 = plot(pltx, plty[:,3])
    display(p3)

end