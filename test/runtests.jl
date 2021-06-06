using Kinetic

# initialization
set, ctr, xface, yface, t = initialize("config.toml")

# solution algorithm
t = solve!(set, ctr, xface, yface, t)

# visualization
plot_contour(set, ctr)
