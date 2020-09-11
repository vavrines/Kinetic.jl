"""
Wave propagation method for the Maxwell's equations

* @param[in]: variables in left-left, left, right, and right-right cells
* @param[in]: eigenmatrix (A), eigenvalue (D)
* @param[in]: speed of light (sol)
* @param[in]: auxiliary parameters (χₑ, νᵦ)

"""
function flux_em!(
    femL::AbstractArray{<:Real,1},
    femR::AbstractArray{<:Real,1},
    ELL::AbstractArray{<:Real,1},
    BLL::AbstractArray{<:Real,1},
    EL::AbstractArray{<:Real,1},
    BL::AbstractArray{<:Real,1},
    ER::AbstractArray{<:Real,1},
    BR::AbstractArray{<:Real,1},
    ERR::AbstractArray{<:Real,1},
    BRR::AbstractArray{<:Real,1},
    ϕL::Real,
    ϕR::Real,
    ψL::Real,
    ψR::Real,
    dxL::Real,
    dxR::Real,
    Ap::AbstractArray{<:Real,2},
    An::AbstractArray{<:Real,2},
    D::AbstractArray{<:Real,1},
    sol::Real,
    χ::Real,
    ν::Real,
    dt::Real,
)

    @assert length(femL) == length(femR) == 8

    slop = zeros(8, 8)
    slop[3, 1] = -0.5 * sol^2 * (BR[2] - BL[2]) + 0.5 * sol * (ER[3] - EL[3])
    slop[5, 1] = 0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 2] = 0.5 * sol^2 * (BR[3] - BL[3]) + 0.5 * sol * (ER[2] - EL[2])
    slop[6, 2] = 0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 3] = 0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 3] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 4] = 0.5 * sol * ν * (BR[1] - BR[1])
    slop[8, 4] = 0.5 * sol^2 * χ * (BR[1] - BR[1])
    slop[3, 5] = -0.5 * sol^2 * (BR[2] - BL[2]) - 0.5 * sol * (ER[3] - EL[3])
    slop[5, 5] = -0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 6] = 0.5 * sol^2 * (BR[3] - BL[3]) - 0.5 * sol * (ER[2] - EL[2])
    slop[6, 6] = -0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 7] = -0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 7] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 8] = -0.5 * sol * ν * (BR[1] - BR[1])
    slop[8, 8] = 0.5 * sol^2 * χ * (BR[1] - BR[1])

    limiter = zeros(8, 8)
    limiter[3, 1] = -0.5 * sol^2 * (BL[2] - BLL[2]) + 0.5 * sol * (EL[3] - ELL[3])
    limiter[5, 1] = 0.5 * sol * (BL[2] - BLL[2]) - 0.5 * (EL[3] - ELL[3])
    limiter[2, 2] = 0.5 * sol^2 * (BL[3] - BLL[3]) + 0.5 * sol * (EL[2] - ELL[2])
    limiter[6, 2] = 0.5 * sol * (BL[3] - BLL[3]) + 0.5 * (EL[2] - ELL[2])
    limiter[1, 3] = 0.5 * sol * χ * (EL[1] - ELL[1])
    limiter[7, 3] = 0.5 * χ * (EL[1] - ELL[1])
    limiter[4, 4] = 0.5 * sol * ν * (BL[1] - BL[1])
    limiter[8, 4] = 0.5 * sol^2 * χ * (BL[1] - BL[1])
    limiter[3, 5] = -0.5 * sol^2 * (BRR[2] - BR[2]) - 0.5 * sol * (ERR[3] - ER[3])
    limiter[5, 5] = -0.5 * sol * (BRR[2] - BR[2]) - 0.5 * (ERR[3] - ER[3])
    limiter[2, 6] = 0.5 * sol^2 * (BRR[3] - BR[3]) - 0.5 * sol * (ERR[2] - ER[2])
    limiter[6, 6] = -0.5 * sol * (BRR[3] - BR[3]) + 0.5 * (ERR[2] - ER[2])
    limiter[1, 7] = -0.5 * sol * χ * (ERR[1] - ER[1])
    limiter[7, 7] = 0.5 * χ * (ERR[1] - ER[1])
    limiter[4, 8] = -0.5 * sol * ν * (BRR[1] - BRR[1])
    limiter[8, 8] = 0.5 * sol^2 * χ * (BRR[1] - BRR[1])

    for i = 1:8
        limiter_theta = sum(slop[:, i] .* limiter[:, i]) / (sum(slop[:, i] .^ 2) + 1.e-7)
        slop[:, i] .*=
            max(0.0, min(min((1.0 + limiter_theta) / 2.0, 2.0), 2.0 * limiter_theta))
    end

    for i = 1:8
        femL[i] =
            sum(An[i, 1:3] .* (ER .- EL)) +
            sum(An[i, 4:6] .* (BR .- BL)) +
            An[i, 7] * (ϕR - ϕL) +
            An[i, 8] * (ψR - ψL) +
            0.5 * sum(
                fortsign.(1.0, D) .* (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D)) .*
                slop[i, :],
            )
        femR[i] =
            sum(Ap[i, 1:3] .* (ER .- EL)) +
            sum(Ap[i, 4:6] .* (BR .- BL)) +
            Ap[i, 7] * (ϕR - ϕL) +
            Ap[i, 8] * (ψR - ψL) -
            0.5 * sum(
                fortsign.(1.0, D) .* (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D)) .*
                slop[i, :],
            )
    end

end


function flux_emx!(
    femL::AbstractArray{<:Real,1},
    femR::AbstractArray{<:Real,1},
    femRU::AbstractArray{<:Real,1},
    femRD::AbstractArray{<:Real,1},
    femLU::AbstractArray{<:Real,1},
    femLD::AbstractArray{<:Real,1},
    ELL::AbstractArray{<:Real,1},
    BLL::AbstractArray{<:Real,1},
    EL::AbstractArray{<:Real,1},
    BL::AbstractArray{<:Real,1},
    ER::AbstractArray{<:Real,1},
    BR::AbstractArray{<:Real,1},
    ERR::AbstractArray{<:Real,1},
    BRR::AbstractArray{<:Real,1},
    ϕL::Real,
    ϕR::Real,
    ψL::Real,
    ψR::Real,
    dxL::Real,
    dxR::Real,
    A1p::AbstractArray{<:Real,2},
    A1n::AbstractArray{<:Real,2},
    A2p::AbstractArray{<:Real,2},
    A2n::AbstractArray{<:Real,2},
    D::AbstractArray{<:Real,1},
    sol::Real,
    χ::Real,
    ν::Real,
    dt::Real,
)

    slop = zeros(8, 8)
    slop[3, 1] = -0.5 * sol^2 * (BR[2] - BL[2]) + 0.5 * sol * (ER[3] - EL[3])
    slop[5, 1] = 0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 2] = 0.5 * sol^2 * (BR[3] - BL[3]) + 0.5 * sol * (ER[2] - EL[2])
    slop[6, 2] = 0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 3] = 0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 3] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 4] = 0.0
    slop[8, 4] = 0.0
    slop[3, 5] = -0.5 * sol^2 * (BR[2] - BL[2]) - 0.5 * sol * (ER[3] - EL[3])
    slop[5, 5] = -0.5 * sol * (BR[2] - BL[2]) - 0.5 * (ER[3] - EL[3])
    slop[2, 6] = 0.5 * sol^2 * (BR[3] - BL[3]) - 0.5 * sol * (ER[2] - EL[2])
    slop[6, 6] = -0.5 * sol * (BR[3] - BL[3]) + 0.5 * (ER[2] - EL[2])
    slop[1, 7] = -0.5 * sol * χ * (ER[1] - EL[1])
    slop[7, 7] = 0.5 * χ * (ER[1] - EL[1])
    slop[4, 8] = 0.0
    slop[8, 8] = 0.0

    limiter = zeros(8, 8)
    limiter[3, 1] = -0.5 * sol^2 * (BL[2] - BLL[2]) + 0.5 * sol * (EL[3] - ELL[3])
    limiter[5, 1] = 0.5 * sol * (BL[2] - BLL[2]) - 0.5 * (EL[3] - ELL[3])
    limiter[2, 2] = 0.5 * sol^2 * (BL[3] - BLL[3]) + 0.5 * sol * (EL[2] - ELL[2])
    limiter[6, 2] = 0.5 * sol * (BL[3] - BLL[3]) + 0.5 * (EL[2] - ELL[2])
    limiter[1, 3] = 0.5 * sol * χ * (EL[1] - ELL[1])
    limiter[7, 3] = 0.5 * χ * (EL[1] - ELL[1])
    limiter[4, 4] = 0.0
    limiter[8, 4] = 0.0
    limiter[3, 5] = -0.5 * sol^2 * (BRR[2] - BR[2]) - 0.5 * sol * (ERR[3] - ER[3])
    limiter[5, 5] = -0.5 * sol * (BRR[2] - BR[2]) - 0.5 * (ERR[3] - ER[3])
    limiter[2, 6] = 0.5 * sol^2 * (BRR[3] - BR[3]) - 0.5 * sol * (ERR[2] - ER[2])
    limiter[6, 6] = -0.5 * sol * (BRR[3] - BR[3]) + 0.5 * (ERR[2] - ER[2])
    limiter[1, 7] = -0.5 * sol * χ * (ERR[1] - ER[1])
    limiter[7, 7] = 0.5 * χ * (ERR[1] - ER[1])
    limiter[4, 8] = 0.0
    limiter[8, 8] = 0.0

    for i = 1:8
        limiter_theta = sum(slop[:, i] .* limiter[:, i]) / (sum(slop[:, i] .^ 2) + 1.e-7)
        slop[:, i] .*=
            max(0.0, min(min((1.0 + limiter_theta) / 2.0, 2.0), 2.0 * limiter_theta))
    end

    for i = 1:8
        femL[i] =
            sum(A1n[i, 1:3] .* (ER .- EL)) +
            sum(A1n[i, 4:6] .* (BR .- BL)) +
            A1n[i, 7] * (ϕR - ϕL) +
            A1n[i, 8] * (ψR - ψL) +
            0.5 * sum(
                fortsign.(1.0, D) .* (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D)) .*
                slop[i, :],
            )
        femR[i] =
            sum(A1p[i, 1:3] .* (ER .- EL)) +
            sum(A1p[i, 4:6] .* (BR .- BL)) +
            A1p[i, 7] * (ϕR - ϕL) +
            A1p[i, 8] * (ψR - ψL) -
            0.5 * sum(
                fortsign.(1.0, D) .* (1.0 .- dt ./ (0.5 * (dxL + dxR)) .* abs.(D)) .*
                slop[i, :],
            )
    end

    # high order correction
    for i = 1:8
        femL[i] +=
            0.5 *
            sum(fortsign.(1.0, D) .* (1.0 .- dt / (dxL + dxR) * abs.(D)) .* slop[i, :])
        femR[i] -=
            0.5 *
            sum(fortsign.(1.0, D) .* (1.0 .- dt / (dxL + dxR) * abs.(D)) .* slop[i, :])
    end

    # transverse correction
    for i = 1:8
        femRU[i] = sum(A2p[i, :] .* femR)
        femRD[i] = sum(A2n[i, :] .* femR)
        femLU[i] = sum(A2p[i, :] .* femL)
        femLD[i] = sum(A2n[i, :] .* femL)
    end
    femRU .*= -0.5 * dt / (dxL + dxR)
    femRD .*= -0.5 * dt / (dxL + dxR)
    femLU .*= -0.5 * dt / (dxL + dxR)
    femLD .*= -0.5 * dt / (dxL + dxR)

end


function flux_emy!(
    femL::AbstractArray{<:Real,1},
    femR::AbstractArray{<:Real,1},
    femRU::AbstractArray{<:Real,1},
    femRD::AbstractArray{<:Real,1},
    femLU::AbstractArray{<:Real,1},
    femLD::AbstractArray{<:Real,1},
    ELL::AbstractArray{<:Real,1},
    BLL::AbstractArray{<:Real,1},
    EL::AbstractArray{<:Real,1},
    BL::AbstractArray{<:Real,1},
    ER::AbstractArray{<:Real,1},
    BR::AbstractArray{<:Real,1},
    ERR::AbstractArray{<:Real,1},
    BRR::AbstractArray{<:Real,1},
    ϕL::Real,
    ϕR::Real,
    ψL::Real,
    ψR::Real,
    dxL::Real,
    dxR::Real,
    A1p::AbstractArray{<:Real,2},
    A1n::AbstractArray{<:Real,2},
    A2p::AbstractArray{<:Real,2},
    A2n::AbstractArray{<:Real,2},
    D::AbstractArray{<:Real,1},
    sol::Real,
    χ::Real,
    ν::Real,
    dt::Real,
)

    slop = zeros(8, 8)
    slop[3, 1] = 0.5 * sol^2 * (BR[1] - BL[1]) + 0.5 * sol * (ER[3] - EL[3])
    slop[4, 1] = 0.5 * sol * (BR[1] - BL[1]) + 0.5 * (ER[3] - EL[3])
    slop[1, 2] = -0.5 * sol^2 * (BR[3] - BL[3]) + 0.5 * sol * (ER[1] - EL[1])
    slop[6, 2] = 0.5 * sol * (BR[3] - BL[3]) - 0.5 * (ER[1] - EL[1])
    slop[2, 3] = 0.5 * sol * χ * (ER[2] - EL[2])
    slop[7, 3] = 0.5 * χ * (ER[2] - EL[2])
    slop[5, 4] = 0.0
    slop[8, 4] = 0.0
    slop[3, 5] = 0.5 * sol^2 * (BR[1] - BL[1]) - 0.5 * sol * (ER[3] - EL[3])
    slop[4, 5] = -0.5 * sol * (BR[1] - BL[1]) + 0.5 * (ER[3] - EL[3])
    slop[1, 6] = -0.5 * sol^2 * (BR[3] - BL[3]) - 0.5 * sol * (ER[1] - EL[1])
    slop[6, 6] = -0.5 * sol * (BR[3] - BL[3]) - 0.5 * (ER[1] - EL[1])
    slop[2, 7] = -0.5 * sol * χ * (ER[2] - EL[2])
    slop[7, 7] = 0.5 * χ * (ER[2] - EL[2])
    slop[5, 8] = 0.0
    slop[8, 8] = 0.0

    limiter = zeros(8, 8)
    limiter[3, 1] = 0.5 * sol^2 * (BL[1] - BLL[1]) + 0.5 * sol * (EL[3] - ELL[3])
    limiter[4, 1] = 0.5 * sol * (BL[1] - BLL[1]) + 0.5 * (EL[3] - ELL[3])
    limiter[1, 2] = -0.5 * sol^2 * (BL[3] - BLL[3]) + 0.5 * sol * (EL[1] - ELL[1])
    limiter[6, 2] = 0.5 * sol * (BL[3] - BLL[3]) - 0.5 * (EL[1] - ELL[1])
    limiter[2, 3] = 0.5 * sol * χ * (EL[2] - ELL[2])
    limiter[7, 3] = 0.5 * χ * (EL[2] - ELL[2])
    limiter[5, 4] = 0.0
    limiter[8, 4] = 0.0
    limiter[3, 5] = 0.5 * sol^2 * (BRR[1] - BR[1]) - 0.5 * sol * (ERR[3] - ER[3])
    limiter[4, 5] = -0.5 * sol * (BRR[1] - BR[1]) + 0.5 * (ERR[3] - ER[3])
    limiter[1, 6] = -0.5 * sol^2 * (BRR[3] - BR[3]) - 0.5 * sol * (ERR[1] - ER[1])
    limiter[6, 6] = -0.5 * sol * (BRR[3] - BR[3]) - 0.5 * (ERR[1] - ER[1])
    limiter[2, 7] = -0.5 * sol * χ * (ERR[2] - ER[2])
    limiter[7, 7] = 0.5 * χ * (ERR[2] - ER[2])
    limiter[5, 8] = 0.0
    limiter[8, 8] = 0.0

    for i = 1:8
        limiter_theta = sum(slop[:, i] .* limiter[:, i]) / (sum(slop[:, i] .^ 2) + 1.e-7)
        slop[:, i] .*=
            max(0.0, min(min((1.0 + limiter_theta) / 2.0, 2.0), 2.0 * limiter_theta))
    end

    for i = 1:8
        femL[i] =
            sum(A2n[i, 1:3] .* (ER .- EL)) +
            sum(A2n[i, 4:6] .* (BR .- BL)) +
            A2n[i, 7] * (ϕR - ϕL) +
            A2n[i, 8] * (ψR - ψL) +
            0.5 * sum(
                fortsign.(1.0, D) .* (1.0 .- dt ./ (0.5 * (dyL + dyR)) .* abs.(D)) .*
                slop[i, :],
            )
        femR[i] =
            sum(A2p[i, 1:3] .* (ER .- EL)) +
            sum(A2p[i, 4:6] .* (BR .- BL)) +
            A2p[i, 7] * (ϕR - ϕL) +
            A2p[i, 8] * (ψR - ψL) -
            0.5 * sum(
                fortsign.(1.0, D) .* (1.0 .- dt ./ (0.5 * (dyL + dyR)) .* abs.(D)) .*
                slop[i, :],
            )
    end

    # high order correction
    for i = 1:8
        femL[i] +=
            0.5 *
            sum(fortsign.(1.0, D) .* (1.0 .- dt / (dyL + dyR) .* abs.(D)) .* slop[i, :])
        femR[i] -=
            0.5 *
            sum(fortsign.(1.0, D) .* (1.0 .- dt / (dyL + dyR) .* abs.(D)) .* slop[i, :])
    end

    # transverse correction
    for i = 1:8
        femRU[i] = sum(A1p[i, :] .* femR)
        femRD[i] = sum(A1n[i, :] .* femR)
        femLU[i] = sum(A1p[i, :] .* femL)
        femLD[i] = sum(A1n[i, :] .* femL)
    end
    femRU .*= -0.5 * dt / (dyL + dyR)
    femRD .*= -0.5 * dt / (dyL + dyR)
    femLU .*= -0.5 * dt / (dyL + dyR)
    femLD .*= -0.5 * dt / (dyL + dyR)

end
