# ------------------------------------------------------------
# Solvers
# ------------------------------------------------------------


export SolverSet


# ------------------------------------------------------------
# Structure of solver setup
# ------------------------------------------------------------
struct SolverSet <: AbstractSolverSet

	# setup
	set :: AbstractSetup

	# physical space
	pSpace :: AbstractPhysicalSpace

	# velocity space
	vSpace :: AbstractVelocitySpace

	# gas property
	gas :: AbstractProperty

	# initial and boundary condition
	ib :: AbstractCondition

	# file system
	outputFolder :: String

	# constructor
	#SolverSet() = SolverSet("./config/config.txt")

	function SolverSet(configfilename::String)
		
		#--- with allowrance ---#
		# read following data from text file
		#=
		allowed = [ "case", "space", "nSpecies", "interpOrder", "limiter", "cfl", "maxTime", 
				    "x0", "x1", "nx", "nxg", "pMeshType",
				    "u0", "u1", "nu", "nug", "vMeshType",
				    "knudsen", "mach", "prandtl", "inK", "omega", "alphaRef", "omegaRef" ]
		D = read_dict(configfilename, allowed)

		# define variables
		case = D["case"]
        space = D["space"]
		interpOrder = D["interpOrder"]
		limiter = D["limiter"]
		cfl = D["cfl"]
		maxTime = D["maxTime"]

		x0 = D["x0"]
		x1 = D["x1"]
		nx = D["nx"]
        nxg = D["nxg"]
        pMeshType = D["pMeshType"]

		u0 = D["u0"]
		u1 = D["u1"]
        nu = D["nu"]
        nug = D["nug"]
		vMeshType = D["vMeshType"]

		knudsen = D["knudsen"]
		mach = D["mach"]
		prandtl = D["prandtl"]
        inK = D["inK"]
        omega = D["omega"]
		alphaRef = D["alphaRef"]
        omegaRef = D["omegaRef"]
		=#

		#--- without allowrance ---#
		D = read_dict(configfilename)

		# automatically generate variables from dictionary
		for key in keys(D)
			s = Symbol(key)
			@eval $s = $(D[key])
		end

		# generate data structure
		dim = parse(Int, space[1])
		gasD = ifelse( parse(Int, space[3]) >= 3, 3, dim ) # in case of plasma
		γ = heat_capacity_ratio(inK, gasD)		
		set = Setup(case, space, nSpecies, interpOrder, limiter, cfl, maxTime)
		
		if dim == 1
			pSpace = PSpace1D(x0, x1, nx, pMeshType, nxg)
		elseif dim == 2
			pSpace = PSpace2D(x0, x1, nx, y0, y1, ny, pMeshType, nxg, nyg)
		else
		end

		if case == "shock"
			μᵣ = ref_vhs_vis(knudsen, alphaRef, omegaRef)
			gas = GasProperty(knudsen, mach, prandtl, inK, γ, omega, alphaRef, omegaRef, μᵣ)

			if space == "1d1f1v"
				vSpace = VSpace1D(u0, u1, nu, vMeshType, nug)

				wL, primL, hL, bcL, wR, primR, hR, bcR = ib_rh(mach, γ, vSpace.u)
				ib = IB1D1F(wL, primL, hL, bcL, wR, primR, hR, bcR)
			elseif space == "1d2f1v"
				vSpace = VSpace1D(u0, u1, nu, vMeshType, nug)

				wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR = ib_rh(mach, γ, vSpace.u, inK)
				ib = IB1D2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
			elseif space == "1d1f3v"
				vSpace = VSpace3D(u0, u1, nu, v0, v1, nv, w0, w1, nw, vMeshType, nug, nvg, nwg)

				wL, primL, fL, bcL, wR, primR, fR, bcR = ib_rh(mach, γ, vSpace.u, vSpace.v, vSpace.w)
				ib = IB1D1F(wL, primL, fL, bcL, wR, primR, fR, bcR)
			end
			
			if space == "1d1f1v"
				wL, primL, hL, bcL,
				wR, primR, hR, bcR = ib_rh(mach, γ, vSpace.u)

				ib = IB1D1F(wL, primL, hL, bcL, wR, primR, hR, bcR)
			elseif space == "1d2f1v"
				wL, primL, hL, bL, bcL,
				wR, primR, hR, bR, bcR = ib_rh(mach, γ, vSpace.u, inK)

				ib = IB1D2F(wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR)
			end
		elseif case == "brio-wu"
			v0 = u0 * sqrt(mi / me)
			v1 = u1 * sqrt(mi / me)
			kne = knudsen * (me / mi)

			vSpace = MVSpace1D(u0, u1, v0, v1, nu, vMeshType, nug)
			gas = PlasmaProperty([knudsen, kne], mach, prandtl, inK, γ, mi, ni, me, ne, lD, rL, sol, echi, bnu)
		
			wL, primL, h0L, h1L, h2L, h3L, bcL, EL, BL, lorenzL, 
            wR, primR, h0R, h1R, h2R, h3R, bcR, ER, BR, lorenzR = ib_briowu(γ, vSpace.u, mi, me)
		    ib = MIB1D4F( wL, primL, h0L, h1L, h2L, h3L, bcL, EL, BL, lorenzL, 
           				  wR, primR, h0R, h1R, h2R, h3R, bcR, ER, BR, lorenzR )
		end

		# create working directory
		identifier = string(Dates.now(), "/")
		#outputFolder = string("../out/", replace(identifier, ":"=>"."))
		outputFolder = replace(identifier, ":"=>".")
		mkdir(outputFolder)
		mkdir(string(outputFolder, "data/"))
		cp(configfilename, string(outputFolder, "config.txt"))

		# create new struct
		new(set, pSpace, vSpace, gas, ib, outputFolder)

	end # function

end # struct


# ------------------------------------------------------------
# Solution algorithm
# ------------------------------------------------------------
function solve!( KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1}, 
				 face::Array{<:AbstractInterface1D,1}, simTime::Float64 )
	
	#--- setup ---#
	dim = ifelse( parse(Int, KS.set.space[3]) >= 3, 3, parse(Int, KS.set.space[1]) ) # dimension

	iter = 0
	dt = 0.
	res = zeros(dim+2)
	#write_jld(KS, ctr, simTime)

	#--- main loop ---#
	#while KS.simTime < KS.maxTime
	while true

		dt = timestep(KS, ctr, simTime)

		reconstruct!(KS, ctr)

		evolve!(KS, ctr, face, dt)

		update!(KS, ctr, face, dt, res)

		iter += 1
		simTime += dt
		
		if iter%1000 == 0
			println("iter: $(iter), time: $(simTime), dt: $(dt), res: $(res[1:dim+2])")

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

end # function


# ------------------------------------------------------------
# Calculate time step
# ------------------------------------------------------------
function timestep(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1}, simTime::Float64)

    tmax = 0.0

	if KS.set.nSpecies == 1

		Threads.@threads for i=1:KS.pSpace.nx
			@inbounds prim = ctr[i].prim
			sos = sound_speed(prim, KS.gas.γ)
			vmax = max(KS.vSpace.u1, abs(prim[2])) + sos
			@inbounds tmax = max(tmax, vmax / ctr[i].dx)
		end

	elseif KS.set.nSpecies == 2

		Threads.@threads for i=1:KS.pSpace.nx
			@inbounds prim = ctr[i].prim
			sos = sound_speed(prim, KS.gas.γ)
			vmax = max(maximum(KS.quad.u1), maximum(abs.(prim[2,:]))) + sos
			@inbounds tmax = ifelse( KS.set.space == "1d4f", 
									 max(tmax, vmax / ctr[i].dx, KS.gas.sol / ctr[i].dx), 
									 max(tmax, vmax / ctr[i].dx) )
		end

	end

    dt = KS.set.cfl / tmax
    dt = ifelse(dt < (KS.set.maxTime - simTime), dt, KS.set.maxTime - simTime)

    return dt

end


# ------------------------------------------------------------
# Reconstruction
# ------------------------------------------------------------
function reconstruct!(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1})

	if KS.set.interpOrder == 1
		return
	end

	#ctr[1].sw .= reconstruct3(ctr[1].w, ctr[2].w, 0.5*(ctr[1].dx+ctr[2].dx))
	#ctr[KS.nx].sw .= reconstruct3(ctr[KS.nx-1].w, ctr[KS.nx].w, 0.5*(ctr[KS.nx-1].dx+ctr[KS.nx].dx))

	# macroscopic variables
	Threads.@threads for i=2:KS.pSpace.nx-1
		@inbounds ctr[i].sw .= reconstruct3( ctr[i-1].w, ctr[i].w, ctr[i+1].w, 
											 0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
											 KS.set.limiter )									
	end
	
	# particle distribution function
	Threads.@threads for i=2:KS.pSpace.nx-1
		if KS.set.space == "1d1f"
			@inbounds ctr[i].sf .= reconstruct3( ctr[i-1].f, ctr[i].f, ctr[i+1].f, 
												 0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												 KS.set.limiter )
		elseif KS.set.space == "1d2f"
			@inbounds ctr[i].sh .= reconstruct3( ctr[i-1].h, ctr[i].h, ctr[i+1].h, 
												 0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												 KS.set.limiter )
			@inbounds ctr[i].sb .= reconstruct3( ctr[i-1].b, ctr[i].b, ctr[i+1].b, 
												 0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												 KS.set.limiter )	
		elseif KS.set.space == "1d4f"
			@inbounds ctr[i].sh0 .= reconstruct3( ctr[i-1].h0, ctr[i].h0, ctr[i+1].h0, 
												  0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												  KS.set.limiter )		
			@inbounds ctr[i].sh1 .= reconstruct3( ctr[i-1].h1, ctr[i].h1, ctr[i+1].h1, 
												  0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												  KS.set.limiter )	
			@inbounds ctr[i].sh2 .= reconstruct3( ctr[i-1].h2, ctr[i].h2, ctr[i+1].h2, 
												  0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												  KS.set.limiter )	
			@inbounds ctr[i].sh3 .= reconstruct3( ctr[i-1].h3, ctr[i].h3, ctr[i+1].h3, 
												  0.5 * (ctr[i-1].dx + ctr[i].dx), 0.5 * (ctr[i].dx + ctr[i+1].dx),
												  KS.set.limiter )							 								 									
	
		end
	end

end


# ------------------------------------------------------------
# Evolution
# ------------------------------------------------------------
function evolve!(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1}, face::Array{Interface1D1F,1}, dt::Float64)

    #if KS.set.case == "heat"
#		flux_maxwell!(KS.ib.bcL, face[1], ctr[1], 1, dt)
#    end

	if KS.set.space == "1d1f"

		Threads.@threads for i=2:KS.pSpace.nx
		#	@inbounds face[i].fw, face[i].ff = flux_kfvs( ctr[i-1].f .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sf, 
		#												  ctr[i].f .- 0.5 .* ctr[i].dx .* ctr[i].sf, 
		#												  KS.vSpace.u, KS.vSpace.weights, dt, ctr[i-1].sf, ctr[i].sf )

			@inbounds face[i].fw, face[i].ff = 
			flux_kcu( ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw, ctr[i-1].f .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sf, 
			ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw, ctr[i].f .- 0.5 .* ctr[i].dx .* ctr[i].sf,
			KS.vSpace.u, KS.vSpace.weights, KS.gas.K, KS.gas.γ, KS.gas.μᵣ, KS.gas.ω, KS.gas.Pr, dt )
		end

	elseif KS.set.space == "1d2f"

		Threads.@threads for i=2:KS.pSpace.nx
			@inbounds face[i].fw, face[i].fh, face[i].fb = 
			flux_kcu( ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw, 
			ctr[i-1].h .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh, 
			ctr[i-1].b .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sb, 
			ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw, 
			ctr[i].h .- 0.5 .* ctr[i].dx .* ctr[i].sh,
			ctr[i].b .- 0.5 .* ctr[i].dx .* ctr[i].sb,
			KS.vSpace.u, KS.vSpace.weights, KS.gas.K, KS.gas.γ, KS.gas.μᵣ, KS.gas.ω, KS.gas.Pr, dt )
		end

	elseif KS.set.space == "1d4f"

		if KS.set.nSpecies == 2
			Threads.@threads for i=2:KS.pSpace.nx
				@inbounds face[i].fw, face[i].fh0, face[i].fh1, face[i].fh2, face[i].fh3 = 
				flux_kcu( ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw, 
						ctr[i-1].h0 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh0, 
						ctr[i-1].h1 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh1, 
						ctr[i-1].h2 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh2, 
						ctr[i-1].h3 .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sh3, 
						ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw, 
						ctr[i].h0 .- 0.5 .* ctr[i].dx .* ctr[i].sh0,
						ctr[i].h1 .- 0.5 .* ctr[i].dx .* ctr[i].sh1,
						ctr[i].h2 .- 0.5 .* ctr[i].dx .* ctr[i].sh2,
						ctr[i].h3 .- 0.5 .* ctr[i].dx .* ctr[i].sh3,
						KS.vSpace.u, KS.vSpace.weights, KS.gas.K, KS.gas.γ, 
						KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.knudsen[1], dt )
			
				@inbounds face[i].femL, face[i].femR = 
				flux_em( ctr[i-2].E, ctr[i-2].B, ctr[i-1].E, ctr[i-1].B, 
						 ctr[i].E, ctr[i].B, ctr[i+1].E, ctr[i+1].B, 
						 ctr[i-1].ϕ, ctr[i].ϕ, ctr[i-1].ψ, ctr[i].ψ, ctr[i-1].dx, ctr[i].dx,
						 KS.gas.A1p, KS.gas.A1n, KS.gas.D1, KS.gas.sol, KS.gas.χ, KS.gas.ν, dt ) 
			end
		end
	
	end

end


# ------------------------------------------------------------
# Update
# ------------------------------------------------------------
function update!(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1}, face::Array{<:AbstractInterface1D,1}, dt::Float64, residual::Array{Float64,1})

	dim = ifelse( parse(Int, KS.set.space[3]) >= 3, 3, parse(Int, KS.set.space[1]) )
    sumRes = zeros(dim+2)
    sumAvg = zeros(dim+2)

    Threads.@threads for i=2:KS.pSpace.nx-1
		@inbounds step!( face[i].fw, face[i].ff, ctr[i].w, ctr[i].prim, ctr[i].f, 
						 face[i+1].fw, face[i+1].ff, KS.gas.γ, KS.vSpace.u, KS.gas.μᵣ, KS.gas.ω,
						 ctr[i].dx, dt, sumRes, sumAvg )
    end

    #if KS.set.case == "heat"
    #    ctr[KS.pSpace.nx].w = deepcopy(ctr[KS.pSpace.nx-1].w)
    #    ctr[KS.pSpace.nx].prim = deepcopy(ctr[KS.pSpace.nx-1].prim)
    #    ctr[KS.pSpace.nx].h = deepcopy(ctr[KS.pSpace.nx-1].h)
    #end

    for i in eachindex(residual)
    	residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end
    
end


# ------------------------------------------------------------
# Stochastic collocation update
# ------------------------------------------------------------
function step!( fwL::Array{Float64,1}, ffL::AbstractArray{Float64,1}, 
				w::Array{Float64,1}, prim::Array{Float64,1}, f::AbstractArray{Float64,1}, 
				fwR::Array{Float64,1}, ffR::AbstractArray{Float64,1}, 
				γ::Float64, u::AbstractArray{Float64,1}, μᵣ::Float64, ω::Float64,
				dx::Float64, dt::Float64, RES::Array{Float64,1}, AVG::Array{Float64,1} )

	#--- store W^n and calculate H^n,\tau^n ---#
	w_old = deepcopy(w)

	#--- update W^{n+1} ---#
	@. w += (fwL - fwR) / dx
	prim .= conserve_prim(w, γ)

	#--- record residuals ---#
	@. RES += (w - w_old)^2
	@. AVG += abs(w)

	#--- calculate M^{n+1} and tau^{n+1} ---#
	M = maxwellian(u, prim)
	τ = vhs_collision_time(prim, μᵣ, ω)

	#--- update distribution function ---#
	for i in eachindex(u)
		f[i] = (f[i] + (ffL[i] - ffR[i]) / dx + dt / τ * M[i]) / (1.0 + dt / τ)
	end

end

function step!( fwL::Array{Float64,1}, fhL::AbstractArray{Float64,1}, fbL::AbstractArray{Float64,1},
				w::Array{Float64,1}, prim::Array{Float64,1}, h::AbstractArray{Float64,1}, b::AbstractArray{Float64,1}, 
				fwR::Array{Float64,1}, fhR::AbstractArray{Float64,1}, fbR::AbstractArray{Float64,1}, 
				K::Float64, γ::Float64, u::AbstractArray{Float64,1}, μᵣ::Float64, ω::Float64,
				dx::Float64, dt::Float64, RES::Array{Float64,1}, AVG::Array{Float64,1} )

	#--- store W^n and calculate H^n,\tau^n ---#
	w_old = deepcopy(w)

	#--- update W^{n+1} ---#
	@. w += (fwL - fwR) / dx
	prim .= conserve_prim(w, γ)

	#--- record residuals ---#
	@. RES += (w - w_old)^2
	@. AVG += abs(w)

	#--- calculate M^{n+1} and tau^{n+1} ---#
	MH = maxwellian(u, prim)
	MB = MH .* K ./ (2. * prim[end])
	τ = vhs_collision_time(prim, μᵣ, ω)

	#--- update distribution function ---#
	for i in eachindex(u)
		h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * MH[i]) / (1.0 + dt / τ)
		b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + dt / τ * MB[i]) / (1.0 + dt / τ)
	end

end


