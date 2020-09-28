
#example
CrossSectionNodes=[2.5 -3.5;2 -4;0 -4;0 4;-2 4;-2.5 3.5]  #member cross-section nodes, in x-y pairs
L=100  #member length
z=0:0.5:L  #z coordinates along the length

#define member deformed shape
a1=L/100*0
a2=L/1000*0
a3=L/1000*10
u=a1*sin.(π*z/L)
v=a2*sin.(π*z/L)
ϕ=a3*sin.(π*z/L)

#mesh the member deformed shape, for undeformed set u=v=ϕ=0
coordinates, connectivity=Mesh3D(CrossSectionNodes, z, u, v, ϕ)

#plot the mesh with Makie
using Makie
scene = poly(coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:false, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 4)
