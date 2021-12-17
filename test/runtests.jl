using Kinetic, Plots

# initialization
set, ctr, xface, yface, t = initialize("config.toml")

# solution algorithm
t = solve!(set, ctr, xface, yface, t)

# visualization
plot(set, ctr)
