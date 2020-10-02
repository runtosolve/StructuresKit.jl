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

# using Plots
# plot(node[:,2], node[:,3], seriestype=:scatter)

A, xc, yc, Ixc, Iyc, Ixyc, Imax, Imin, Th_p, Cw, J, Xs, Ys, w, Bx, By, B1, B2, Qx, Qy, Ai = CrossSection.CUFSMproperties(node, elem)

Q = Qx[1:5]
dA = Ai[1:4]
V=2.5
I = Ixc


#shear stress in flange
τ = CrossSection.shear_stress.(V, Q, I, t)
#integrate stress over flange
Vfx = CrossSection.shear_force(τ, dA)


#shear stress in web
Q=Qx[5:9]
dA = Ai[5:8]
τweb = CrossSection.shear_stress.(V, Q, I, t)
Vweb = CrossSection.shear_force(τweb, dA)
