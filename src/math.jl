# ============================================================
# Mathematical Methods
# ============================================================


export heaviside,
       fortsign,
       solve_linear


# ------------------------------------------------------------
# Heaviside step function
# ------------------------------------------------------------
heaviside(x::Union{Int, AbstractFloat}) = ifelse(x >= 0, 1., 0.)


# ------------------------------------------------------------
# Fortran sign()
# ------------------------------------------------------------
fortsign(x, y) = abs(x) * sign(y)


# ------------------------------------------------------------
# Alternative solver for Linear system
# ------------------------------------------------------------
function solve_linear(A::Array{Float64,2}, b::Array{Float64,1})

    A_temp = deepcopy(A)
    b_temp = deepcopy(b)
    x = similar(b)

    for i=1:length(b)
        temp = maximum(A_temp[:,i])
        for j=i+1:length(b)
            if A_temp[j,i] == temp
                exchangea = A_temp[j,:]
                A_temp[j,:] = A_temp[i,:]
                A_temp[i,:] = exchangea
                exchangeb = b_temp[j]
                b_temp[j] = b_temp[i]
                b_temp[i] = exchangeb
            end
        end
        for j=i+1:length(b)
            temp = A_temp[j,i] / A_temp[i,i]
            for k=i:length(b)
                A_temp[j,k] = A_temp[j,k] - temp * A_temp[i,k]
            end
            b_temp[j] = b_temp[j] - temp * b_temp[i]
        end
        if abs(A_temp[i,i]) < 1.e-5
            println("warning: linear solver zero denominator")
        end
    end

    for i=length(b):-1:1
        for j=i-1:-1:1
            temp = A_temp[j,i] / A_temp[i,i]
            for k=i:length(b)
                A_temp[j,k] = A_temp[j,k] - temp * A_temp[i,k]
            end
            b_temp[j] = b_temp[j] - temp * b_temp[i]
        end
    end

    for i=1:length(b)
        x[i] = b_temp[i] / A_temp[i,i]
        if isnan(x[i])
            println("error: linear solver NaN solution $(i)")
        end
    end

    return x
    
end