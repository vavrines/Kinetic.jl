"""
Kinetic.jl: A Portable Framework for Scientific and Neural Computing

Copyright (c) 2020-2024 Tianbai Xiao (tianbaixiao@gmail.com)
"""
module Kinetic

export 転

"""
轻量化的科学与神经网络计算框架

"転" means "rolling" in Japanese
"""
const 転 = Kinetic

using Reexport
@reexport using KitBase
@reexport using KitML

end
