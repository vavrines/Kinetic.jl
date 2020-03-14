using Kinetic
#cd("D:\\Github\\Kinetic.jl\\example\\")

ks, ctr, face, simTime = initialize("shock.txt");

solve!(ks, ctr, face, simTime)

plot_line(ks, ctr)