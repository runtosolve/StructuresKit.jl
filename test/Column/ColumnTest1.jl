using StructuresKit

shape_name = "W14X90"
#define discretization along member
zcoords = 0:1:20
n=(4, 4, 4, 4, 4)

shape_info = StructuresKit.CrossSection.AISC(shape_name)

#A Ix Iy J Cw xc yc xs ys
section_properties = [(shape_info.A, shape_info.Ix, shape_info.Iy, shape_info.J, shape_info.Cw, shape_info.d/2, 0.0, shape_info.d/2, 0.0)]

#E  ν
material_properties = [(29000,0.30)]

#kx ky kϕ hx hy
springs = [0.0, 0.0, 0.0, 0.0, 0.0]

#member information
#L dL SectionProperties MaterialProperties Springs
member_definitions = [(13.0*12, (13*12)/12,1,1,1)]

num_nodes = Int(member_definitions[1][1]/member_definitions[1][2]) + 1

#P 
loads = 4000.0 * ones(num_nodes)

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 13.0*12]

#imperfection
L = member_definitions[1][1]
dL = member_definitions[1][2]
z = 0: dL :L
Δo = L / 1000 * sin.(π*z/L) 

uo = Δo  #in the weak axis
vo = zeros(Float64, num_nodes)
ϕo = zeros(Float64, num_nodes)

imperfections = [uo, vo, ϕo]

P = 0:100:4000

u_total_predicted_max = zeros(Float64, length(P))
u_total_theor_max = zeros(Float64, length(P))

#Calculate weak axis flexural buckling load
E = material_properties[1][1]
Iy = shape_info.Iy
Pcry = π^2*E*Iy/L^2   #weak axis buckling loads


for i=1:length(P)

    loads = P[i] * ones(num_nodes)
    u, v, ϕ, properties = Column.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

    u_total_predicted_max[i] = maximum(u) + L/1000

    #displacement at midheight
    #Chajes Eq. 1.70
    u_total_theor_max[i] = L/1000*(1/(1-P[i]/Pcry))

end

# using Plots
# plot(u_total_predicted_max, P, markershape = :o)
# plot!(u_total_theor_max, P, markershape = :hexagon)


#calculate percent error at every load compared to Mathematica solution
u_difference= abs.((u_total_predicted_max .- u_total_theor_max))./ u_total_theor_max

MaxError=maximum(u_difference[2:end])
#accept 7% error from numerical solution
#solution becomes less consistent with Chajes as P approaches Pcry
#this could be because the Chajes solution is an approximation, hard to tell if computed solution or Chajes is more correct

@test MaxError <= 0.07


