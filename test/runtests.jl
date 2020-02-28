using Kinetic

prim = [1., 0., 1.]
inK = 0
γ = heat_capacity_ratio(inK, 1) 

w0 = prim_conserve(prim, γ)
prim0 = conserve_prim(w0, γ)
