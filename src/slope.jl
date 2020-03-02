# ------------------------------------------------------------
# Interpolation and Slope Limiters
# ------------------------------------------------------------


export vanleer,
       minmod


vanleer(sL, sR) = 
(fortsign(1., sL) + fortsign(1., sR)) * abs(sL) * abs(sR) / (abs(sL) + abs(sR) + 1.e-7)


minmod(sL, sR) = 
0.5 * (fortsign(1., sL) + fortsign(1., sR)) * min(abs(sR), abs(sL))