using Kinetic

cd(@__DIR__)
ks, ctr, face, t = Kinetic.initialize("sod.txt")

t = Kinetic.solve!(ks, ctr, face, t)

# it's equivalent as the following process
dt = Kinetic.timestep(ks, ctr, t)
nt = Int(floor(ks.set.maxTime / dt))
res = zeros(3)
for iter = 1:nt
    Kinetic.reconstruct!(ks, ctr)
    Kinetic.evolve!(ks, ctr, face, dt)
    Kinetic.update!(ks, ctr, face, dt, res)
end

Kinetic.plot_line(ks, ctr)
