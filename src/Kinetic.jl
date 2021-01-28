"""
Kinetic.jl: A Portable Framework for Computational Fluid Dynamics and Scientific Machine Learning

Copyright (c) 2021 Tianbai Xiao (tianbaixiao@gmail.com)
"""
module Kinetic

export 転

"""
Portable Framework for Computational Fluid Dynamics and Scientific Machine Learning

轻量化的计算流体力学建模和计算框架

"転" means "rolling" in Japannese
"""
const 転 = Kinetic

using Reexport
@reexport using KitBase
@reexport using KitML

end
