# ------------------------------------------------------------
# Solvers
# ------------------------------------------------------------


export SolverSet1D,
	   solve!,
	   timestep,
	   reconstruct!,
	   evolve!,
	   update!,
	   step!


# ------------------------------------------------------------
# Structure of solver setup
# ------------------------------------------------------------
struct SolverSet1D <: AbstractSolverSet

	# setup
	set :: AbstractSetup

	# physical space
	pMesh :: AbstractPhysicalMesh

	# velocity space
	uMesh :: AbstractVelocityMesh

	# gas property
	gas :: AbstractProperty

	# initial and boundary condition
	ib :: AbstractCondition

	# file system
	outputFolder :: String

	#--- initialization ---#
	function SolverSet1D(configfilename::String)
		
		# read following data from text file
		allowed = [ "case", "space", "interpOrder", "limiter", "cfl", "maxTime", 
				    "x0", "x1", "nx", "nxg", "pMeshType",
				    "u0", "u1", "nu", "nug", "vMeshType",
				    "knudsen", "mach", "prandtl", "inK", "omega", "alphaRef", "omegaRef" ]
		D = read_dict(configfilename, allowed)

		# read configuration from text file
		case = D["case"]
        space = D["space"]
		interpOrder = parse(Int64, D["interpOrder"])
		limiter = D["limiter"]
		cfl = parse(Float64, D["cfl"])
		maxTime = parse(Float64, D["maxTime"])

		x0 = parse(Float64, D["x0"])
		x1 = parse(Float64, D["x1"])
		nx = parse(Int64, D["nx"])
        nxg = parse(Int64, D["nxg"])
        pMeshType = D["pMeshType"]

		u0 = parse(Float64, D["u0"])
		u1 = parse(Float64, D["u1"])
        nu = parse(Int64, D["nu"])
        nug = parse(Int64, D["nug"])
		vMeshType = D["vMeshType"]

		Kn = parse(Float64, D["knudsen"])
		Ma = parse(Float64, D["mach"])
		Pr = parse(Float64, D["prandtl"])
        K = parse(Float64, D["inK"])
        ω = parse(Float64, D["omega"])
		αᵣ = parse(Float64, D["alphaRef"])
        ωᵣ = parse(Float64, D["omegaRef"])

		# deduce configuration from existed data
		γ = heat_capacity_ratio(K, parse(Int,space[1]))
        μᵣ = ref_vhs_vis(Kn, αᵣ, ωᵣ)

		# generate data structure
		set = Setup(case, space, interpOrder, limiter, cfl, maxTime)
		pMesh = PMesh1D(x0, x1, nx, pMeshType, nxg)
        uMesh = VMesh1D(u0, u1, nu, vMeshType, nug)
	    gas = GasProperty(Kn, Ma, Pr, K, γ, ω, αᵣ, ωᵣ, μᵣ)

		if case == "shock"
			if space == "1d1f"
				wL, primL, hL, bcL,
				wR, primR, hR, bcR = ib_rh(Ma, γ, uMesh.u)

				ib = IB1D1F(wL, primL, hL, bcL, wR, primR, hR, bcR)
			elseif space == "1d2f"
				wL, primL, hL, bL, bcL,
				wR, primR, hR, bR, bcR = ib_rh(Ma, γ, uMesh.u, K)

				ib = IB1D2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
			end
		else
		end

		# create working directory
		identifier = string(Dates.now(), "/")
		#outputFolder = string("../out/", replace(identifier, ":"=>"."))
		outputFolder = replace(identifier, ":"=>".")
		mkdir(outputFolder)
		mkdir(string(outputFolder, "data/"))
		cp(configfilename, string(outputFolder, "config.txt"))

		# create new struct
		new(set, pMesh, uMesh, gas, ib, outputFolder)

	end # function

end # struct


# ------------------------------------------------------------
# Solution algorithm
# ------------------------------------------------------------
function solve!(KS::SolverSet1D, ctr, face, simTime)
	
	#--- setup ---#
	dim = parse(Int, KS.set.space[1])

	iter = 0
	dt = 0.
	res = zeros(dim+2)
	
	#--- main loop ---#
	#while KS.simTime < KS.maxTime
#=	while true

		dt = timestep(KS, ctr, simTime)

		reconstruct!(KS, ctr)

		evolve!(KS, ctr, face, dt)

		update!(KS, ctr, face, dt, res, i)

		iter += 1
		simTime += dt
		
		if iter%1000 == 0
			println("iter: $(iter), time: $(simTime), dt: $(dt), res: $(res[1:KS.set.dim+2])")

			#if iter%400 == 0
				#write_jld(KS, ctr, iter)
			#end
		end

		if maximum(res) < 5.e-7 || simTime > KS.set.maxTime
		#if simTime > KS.set.maxTime
			break
		end

	end # while loop

	write_jld(KS, ctr, simTime)

	return ctr
=#
end # function


# ------------------------------------------------------------
# Calculate time step
# ------------------------------------------------------------
function timestep(KS::SolverSet1D, ctr::AbstractArray{<:AbstractControlVolume1D,1}, simTime::Float64)

    tmax = 0.0

    Threads.@threads for i=1:KS.mesh.nx
        @inbounds prim = ctr[i].prim
        sos = sound_speed(prim, KS.gas.γ)
        vmax = max(KS.quad.u1, abs(prim[2])) + sos
        @inbounds tmax = max(tmax, vmax / ctr[i].dx)
    end

    dt = KS.set.cfl / tmax
    dt = ifelse(dt < (KS.set.maxTime - simTime), dt, KS.set.maxTime - simTime)

    return dt

end


# ------------------------------------------------------------
# Reconstruction
# ------------------------------------------------------------
function reconstruct!(KS::SolverSet1D, ctr::AbstractArray{<:AbstractControlVolume1D,1})

	if KS.set.interpOrder == 1
		return
	elseif KS.set.interpOrder == 2
		#ctr[1].sw .= reconstruct3(ctr[1].w, ctr[2].w, 0.5*(ctr[1].dx+ctr[2].dx))
		#ctr[KS.nx].sw .= reconstruct3(ctr[KS.nx-1].w, ctr[KS.nx].w, 0.5*(ctr[KS.nx-1].dx+ctr[KS.nx].dx))

	    Threads.@threads for i=2:KS.mesh.nx-1
			@inbounds ctr[i].sw .= reconstruct3( ctr[i-1].w, ctr[i].w, ctr[i+1].w, 
												 0.5*(ctr[i-1].dx+ctr[i].dx), 0.5*(ctr[i].dx+ctr[i+1].dx),
									   			 KS.set.limiter )
	    end
	end

end


# ------------------------------------------------------------
# Evolution
# ------------------------------------------------------------
function evolve!(KS::SolverSet1D, ctr::AbstractArray{<:AbstractControlVolume1D,1}, face::Array{<:AbstractInterface1D,1}, dt::Float64)

    #if KS.set.case == "heat"
#		flux_maxwell!(KS.ib.bcL, face[1], ctr[1], 1, dt)
#    end

    Threads.@threads for i=2:KS.mesh.nx
		@inbounds flux_kfvs( ctr[i-1].h, ctr[i].h, KS.vMesh.u, KS.vMesh.weights, dt, ctr[i-1].sh, ctr[i].h )
    end

end


# ------------------------------------------------------------
# Update
# ------------------------------------------------------------
function update!(KS::SolverSet1D, ctr::AbstractArray{<:AbstractControlVolume1D,1}, face::Array{<:AbstractInterface1D,1}, dt::AbstractFloat, residual::Array{Float64,1})

	dim = parse(Int, KS.set.space[1])
    sumRes = zeros(dim+2)
    sumAvg = zeros(dim+2)

    Threads.@threads for i=2:KS.mesh.nx-1
		@inbounds step!( face[i].fw, face[i].fh, ctr[i].w, ctr[i].prim, ctr[i].h, ctr[i].dx,
						 face[i+1].fw, face[i+1].fh, KS.gas.γ, KS.vMesh.u, KS.gas.μᵣ, KS.gas.ω,
						 dt, sumRes, sumAvg )
    end

    #if KS.set.case == "heat"
    #    ctr[KS.mesh.nx].w = deepcopy(ctr[KS.mesh.nx-1].w)
    #    ctr[KS.mesh.nx].prim = deepcopy(ctr[KS.mesh.nx-1].prim)
    #    ctr[KS.mesh.nx].h = deepcopy(ctr[KS.mesh.nx-1].h)
    #end

    for i in eachindex(residual)
    	residual[i] = sqrt(sumRes[i] * KS.mesh.nx) / (sumAvg[i] + 1.e-7)
    end
    
end


# ------------------------------------------------------------
# Stochastic collocation update
# ------------------------------------------------------------
function step!( fwL::Array{Float64,1}, fhL::AbstractArray{Float64,1}, 
				w::Array{Float64,1}, prim::Array{Float64,1}, h::AbstractArray{Float64,1}, dx::Float64, 
				fwR::Array{Float64,1}, fhR::AbstractArray{Float64,1}, 
				γ::Float64, u::AbstractArray{Float64,1}, μᵣ::Float64, ω::Float64,
				dt::Float64, RES::Array{Float64,1}, AVG::Array{Float64,1} )

	#--- store W^n and calculate H^n,\tau^n ---#
	w_old = deepcopy(w)

	#--- update W^{n+1} and calculate H^{n+1}, B^{n+1}, \tau^{n+1} ---#
	@. w += (fwL - fwR) / dx
	prim = conserve_prim(w, γ)

	#--- record residuals ---#
	@. RES += (w - w_old)^2
	@. AVG += abs(w)

	#--- construct point collocations of H^{n+1} ---#
	H = maxwellian(u, prim)

	τ = vhs_collision_time(prim, μᵣ, ω)

	#--- update distribution function ---#
	for i in eachindex(u)
		h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * H[i]) / (1.0 + dt / τ)
	end

end