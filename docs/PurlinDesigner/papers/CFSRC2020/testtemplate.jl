using StructuresKit

#bring in cross-section dimensions
#Z200-D2

Bc = 68.4
Dc = 26.4
αc = 92
θc = 47

Bt = 72.5
Dt = 27.5
α = 92
θt = 53

H = 202
r = 6.6
t = 2.59

Fy = 420  #MPa

#calculate centerline dimensions

CorZ = 2    #Z section considered here
center = 0  #outside dimensions
kipin = 0

nh = 4
nb1 = 4
nb2 = 4
nd1 = 4
nd2 = 4
nr1 = 4
nr2 = 4
nr3 = 4
nr4 = 4


prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,Bt,Bc,Dt,Dc,r,r,r,r,θt,θc,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

using Plots
plot(node[:,2], node[:, 3], markershape = :circle)


kϕ = 1339  #N-mm/rad/mm   c=34 mm closest to test 
