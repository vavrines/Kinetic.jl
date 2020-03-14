using Kinetic: AbstractPhysicalMesh, flux_kfvs, PMesh1D, heaviside, fortsign, vanleer, minmod,
               gauss_moments, moments_conserve, moments_conserve_slope, discrete_moments, maxwellian, 
               conserve_prim, prim_conserve, heat_capacity_ratio, ref_vhs_vis, sound_speed, vhs_collision_time, 
               aap_hs_collision_time, aap_hs_prim, newton_cotes

#cd("D:\\Github\\Kinetic.jl\\example\\")

prim = [1, 0, 1]
w = prim_conserve(prim, 1.4)