using StructuresKit


#Confirm that VX is the first moment area about the x axis in the CrossSection.CUFSMProperties function.

#Use Beer and Johnston Mechanics of Materials, Example 6.06

#Cee section without lips

CorZ = 1
t = 0.15
H = 6.0
B = 4.0
D = 0.0
θ = 0.0


r = 0.0
kipin = 0
center = 1
n = 4

prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,B,B,D,D,r,r,r,r,θ,θ,t,n,n,n,n,n,n,n,n,n,kipin,center)

using Plots
plot(node[:,2], node[:,3], seriestype=:scatter)

A, xc, yc, Ixc, Iyc, Ixyc, Imax, Imin, Th_p, Cw, J, Xs, Ys, w, Bx, By, B1, B2, VX, VY = CrossSection.CUFSMproperties(node, elem)


#Stress at bottom flange-web intersection is 2.22 ksi (from left to right) for V=2.5 kips down.

#VX is the same as the first moment area taken about the x axis, call it Qx

τ_solution = 2.22

Q = VX[5]
V = -2.5
τ =  V*Q/(Ixc*t)
