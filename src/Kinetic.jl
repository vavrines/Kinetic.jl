"""
Kinetic.jl: A Lightweight Framework for Kinetic Theory and Scientific Machine Learning

Copyright (c) 2020 Tianbai Xiao <<tianbaixiao@gmail.com>>
"""
module Kinetic

export 転

"""
Lightweight Framework for Kinetic Theory and Scientific Machine Learning

轻量化的动理论建模和计算框架

"転" means "rolling" in Japannese
"""
const 転 = Kinetic

using Reexport
@reexport using KitBase
@reexport using KitML

end
