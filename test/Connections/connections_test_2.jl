using StructuresKit

#User inputs
t1 = 0.78
t2 = 2.54
fy1 = 150.4
fy2 = 422.4
fu1 = 312.1
fu2 = 534.2
screw_diameter = 4.75

D1, D2, D3, D4, P1, P2, P3, P4 = Connections.cfs_load_deformation_interpolation_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)

using Plots
plot([0.0, D1, D2, D3, D4], [0.0, P1, P2, P3, P4], markershape = :o)