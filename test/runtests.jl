using Kinetic, Plots

# initialization
set, ctr, xface, yface, t = Kinetic.initialize("config.toml")

# solution algorithm
t = Kinetic.solve!(set, ctr, xface, yface, t)

# visualization
plot(set, ctr)

# short name
転 == Kinetic
