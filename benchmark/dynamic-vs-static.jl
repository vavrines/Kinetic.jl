"""
This file compares the performance by using dynamic and static arrays to store solutions.
"""

using Kinetic, StaticArrays, BenchmarkTools, OffsetArrays
cd(@__DIR__)

#--- dynamic ---#
ks, ctr, face, t = Kinetic.initialize("bench_ds.txt")

@btime Kinetic.solve!(ks, ctr, face, t)
"""The result on my NUC is around 707 ms."""

#--- static ---#
ctr2 = deepcopy(ctr)
for i in eachindex(ctr2)
    if i <= ks.pSpace.nx ÷ 2
        ctr2[i] = ControlVolume1D2F(
            ks.pSpace.x[i],
            ks.pSpace.dx[i],
            MVector{length(ks.ib.wL)}(ks.ib.wL),
            MVector{length(ks.ib.wL)}(ks.ib.primL),
            MVector{length(ks.ib.hL)}(ks.ib.hL),
            MVector{length(ks.ib.hL)}(ks.ib.bL),
        )
    else
        ctr2[i] = ControlVolume1D2F(
            ks.pSpace.x[i],
            ks.pSpace.dx[i],
            MVector{length(ks.ib.wL)}(ks.ib.wR),
            MVector{length(ks.ib.wL)}(ks.ib.primR),
            MVector{length(ks.ib.hL)}(ks.ib.hR),
            MVector{length(ks.ib.hL)}(ks.ib.bR),
        )
    end
end

face2 = deepcopy(face)
for i = 1:ks.pSpace.nx+1
    face2[i] = Interface1D2F(
        MVector{length(ks.ib.wL)}(ks.ib.wL),
        MVector{length(ks.ib.hL)}(ks.ib.hL),
    )
end

@btime Kinetic.solve!(ks, ctr2, face2, t)
"""
The result on my NUC is around 523 ms.
Further improvements are possible if we make all structs (ks, ctr, face) static.
"""
