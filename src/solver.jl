# ------------------------------------------------------------
# Solvers
# ------------------------------------------------------------


export SolverSet1D


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
		allowed = [ "isNewStart", "case", "space", "interpOrder", "cfl", "maxTime", 
				    "x0", "x1", "nx", "nxg", "pMeshType",
				    "u0", "u1", "nu", "nug", "vMeshType",
				    "knudsen", "mach", "prandtl", "inK", "omega", "alphaRef", "omegaRef" ]
		D = readtodict(configfilename, allowed)

		# read configuration from text file
        isNewStart = Bool(parse(Int64, D["isNewStart"]))
		case = D["case"]
        space = D["space"]
        interpOrder = parse(Int64, D["interpOrder"])
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
		set = Setup(isNewStart, case, space, interpOrder, cfl, maxTime)
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
		outputFolder = string("../out/", replace(identifier, ":"=>"."))
		mkdir(outputFolder)
		mkdir(string(outputFolder, "data/"))
		cp(configfilename, string(outputFolder, "config.txt"))

		# create new struct
		new(set, pMesh, uMesh, gas, ib, outputFolder)

	end # function

end # struct