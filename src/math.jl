# ============================================================
# Mathematical Operations
# ============================================================


export global_frame, 
       local_frame


function global_frame(w::Array{Float64,1}, cosa::AbstractFloat, sina::AbstractFloat) 
    
    if length(w) == 2
        G = [ w[1] * cosa - w[2] * sina, w[1] * sina + w[2] * cosa ]
    elseif length(w) == 4
        G = [ w[1], w[2] * cosa - w[3] * sina, w[2] * sina + w[3] * cosa, w[4] ]
    else
        println("error: local -> global")
    end

    return G

end


function local_frame(w::Array{Float64,1}, cosa::AbstractFloat, sina::AbstractFloat)
    
    if length(w) == 2
        L = [ w[1] * cosa + w[2] * sina, w[2] * cosa - w[1] * sina]
    elseif length(w) == 4
        L = [ w[1], w[2] * cosa + w[3] * sina, w[3] * cosa - w[2] * sina, w[4] ]
    else
        println("error: global -> local")
    end

    return L

end