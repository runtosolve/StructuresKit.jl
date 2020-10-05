
using Statistics
using LinearAlgebra

coord = [0.0 0.0
         0.34 0.76
         0.98 0.23]

ends = [1 2 0.01
        2 3 0.01]

A,xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,Cw,B1,B2,wn = cutwp_prop2(coord,ends)
