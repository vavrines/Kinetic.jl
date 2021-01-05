"""
Kinetic.jl: A Lightweight Framework for Computational Fluid Dynamics and Scientific Machine Learning

Copyright (c) 2021 Tianbai Xiao (tianbaixiao@gmail.com)
"""
module Kinetic

if VERSION < v"1.3"
    @warn "Kinetic.jl matches perfectly with Julia 1.3 and newer versions."
end

export 転

"""
Lightweight Framework for Computational Fluid Dynamics and Scientific Machine Learning

轻量化的计算流体力学建模和计算框架

"転" means "rolling" in Japannese
"""
const 転 = Kinetic

using Reexport
@reexport using KitBase
@reexport using KitML

end
