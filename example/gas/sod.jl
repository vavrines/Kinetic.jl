using KitBase, Plots

cd(@__DIR__)
ks, ctr, face, t = initialize("sod.txt")
t = solve!(ks, ctr, face, t)

plot(ks, ctr)
