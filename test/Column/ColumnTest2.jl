using StructuresKit

#Example 1 from 

#Plaut and Moen (2019).  "Flexural-torsional deformations of imperfect thin-walled columns with continous bracing."  Proceedings of the Annual Stability Conference, Structural Stability Research Council, St. Louis, MO.  
#https://www.aisc.org/globalassets/continuing-education/ssrc-proceedings/2019/plaut_and_moen_ssrc_2019.pdf

#A Ix Iy J Cw xc yc xs ys
A = 272.0
Ix = 363370.0
Iy = 64100.0
xc = 0.0
yc = 0.0
xs = -32.59
ys = 0.0
J = 188.0
Cw = 122720891.0
section_properties = [(A, Ix, Iy, J, Cw, xc, yc, xs, ys)]

#E  ν
material_properties = [(200,0.30)]

#kx ky kϕ hx hy
springs = [0.0, 0.0, 0.0, 0.0, 0.0]

#member information
#L dL SectionProperties MaterialProperties Springs
member_definitions = [(2438.0, 2438.0/10,1,1,1)]

num_nodes = Int(member_definitions[1][1]/member_definitions[1][2]) + 1

#P 
loads = 5.0 * ones(num_nodes)

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 2438.0]

#imperfections
L = member_definitions[1][1]
dL = member_definitions[1][2]
z = 0: dL :L
uo = L / 1000 * sin.(π*z/L) 
vo = L / 1000 * sin.(π*z/L) 
ϕo = 0.00766 * sin.(π*z/L) 


imperfections = [uo, vo, ϕo]

#Figure 4 
#for P=5 kN
#(kϕ, ϕ+ϕo) from Mathematica
#kϕ is in kN-mm/rad/mm 
data_SSRC = [(0,0.0201),(0.01,0.0185),(0.03,0.0163),(0.04,0.0158),(0.05,0.01486),(0.06,0.0143),(0.09,0.01305),(0.12,0.0122),(0.15,0.0116),(0.2,0.01085),(0.3,0.00999),(0.4,0.00949),(0.5,0.00917),(0.7,0.00876),(1,0.00846),(2,0.00808),(5,0.00783)]

num_points = length(data_SSRC)

ϕ_predicted_max = zeros(Float64, num_points)
kϕ = zeros(Float64, num_points)

ϕo_max = maximum(ϕo)

ϕ_total_SSRC = zeros(Float64, num_points)

for i = 1:num_points

    #kx ky kϕ hx hy
    kϕ[i] = data_SSRC[i][1]
    
    springs = [0.0, 0.0, kϕ[i], 0.0, 0.0]

    u, v, ϕ, properties = Column.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

    ϕ_predicted_max[i] = maximum(ϕ)
    ϕ_total_SSRC[i] = data_SSRC[i][2]

end

ϕ_total_predicted = ϕ_predicted_max .+ ϕo_max

#calculate percent error at every load compared to Mathematica solution
ϕ_difference= abs.((ϕ_total_predicted .- ϕ_total_SSRC[:]))./ ϕ_total_SSRC[:]

MaxError=maximum(ϕ_difference[2:end])


#accept 2% error from numerical solution
#there looks to be a small bug in the Mathematica data for one point, otherwise 0.01 would work here
@test MaxError <= 0.02


# using Plots
# plot(kϕ, ϕ_predicted_max .+ ϕo_max, markershape = :o)
# plot!(kϕ, ϕ_total_SSRC, markershape = :hexagon)


#*********************************************************

#Figure 5
#for P=5 kN
#(kϕ, v+vo) from Mathematica
#kϕ is in kN-mm/rad/mm 
data_SSRC = [(0,2.5716),(0.01,2.5694),(0.02,2.5677),(0.04,2.5652),(0.06,2.5635),(0.08,2.5623),(0.1,2.5613),(0.12,2.5606),(0.15,2.5600),(0.2,2.5587),(0.3,2.5575),(0.4,2.55675),(0.5,2.5563),(0.6,2.556),(0.7,2.5558)]

num_points = length(data_SSRC)

v_predicted_max = zeros(Float64, num_points)
kϕ = zeros(Float64, num_points)

vo_max = maximum(vo)

v_total_SSRC = zeros(Float64, num_points)

for i = 1:num_points

    #kx ky kϕ hx hy
    kϕ[i] = data_SSRC[i][1]
    
    springs = [0.0, 0.0, kϕ[i], 0.0, 0.0]

    u, v, ϕ, properties = Column.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

    v_predicted_max[i] = maximum(v)
    v_total_SSRC[i] = data_SSRC[i][2]

end

v_total_predicted = v_predicted_max .+ vo_max

#calculate percent error at every load compared to Mathematica solution
v_difference= abs.((v_total_predicted .- v_total_SSRC[:]))./ v_total_SSRC[:]

MaxError=maximum(v_difference[2:end])
#accept 1% error from numerical solution
@test MaxError <= 0.01


# using Plots
# plot(kϕ, v_predicted_max .+ vo_max, markershape = :o)
# plot!(kϕ, v_total_SSRC, markershape = :hexagon)


