# ============================================================
# Methods of Input
# ============================================================


export readtodict


# ------------------------------------------------------------
# Read text into dictionary
#
# >@param[in]  filename     :  configuration text file
# >@param[in]  allowed      :  keywords
# >@return  vars            :  dictionary with values of variables
# ------------------------------------------------------------
function readtodict(filename::String, allowed)

    println("")
    f = open(filename)
    vars = Dict{String, Any}()
    println("Reading config from $filename")
    println("")

    for line in eachline(f)
        # skip comments
        if length(line) == 0 || line[1] == '#' 
            #println("skip comment line")
            continue
        end

        #print("\t");
        println(line)

        var, val = split(line, "=")
        stripped = strip(var)
        if stripped in allowed
            #vars[stripped] = parse(Float64, val)
            vars[stripped] = strip(val)
        end
    end

    println("")
    return vars

end