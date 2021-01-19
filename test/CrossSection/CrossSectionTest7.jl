
using StructuresKit
using DiffEqOperators


shape_name = "W14X90"

shape_info = StructuresKit.CrossSection.AISC(shape_name)

#A Ix Iy J Cw xc yc xs ys
section_properties = [(shape_info.A, shape_info.Ix, shape_info.Iy, shape_info.J, shape_info.Cw, 0.0, shape_info.d/2, 0.0, shape_info.d/2)]

#define column cross-section mesh
n_Wshape=(4, 2, 4, 4, 4)
mesh_size = 0.01
xcoords, ycoords = CrossSection.wshape_nodes(shape_info, n_Wshape)
section_geometry = Column.SectionGeometry(xcoords, ycoords, mesh_size)

#define material properties
Es = 29000.0
σy = 50.0                   #steel yield stress
σy1 = 50.0                    #steel yield stress
σu = 65.3                   #steel ultimate stress
σf= σu      #steel fracture stress#
ϵy = σy / Es              #steel yield strain
ϵy1 = σy / Es * 10    #steel strain at end of yield plateau
ϵu = 0.18             #steel ultimate strain
ϵf = 0.21             #steel fracture strain
n = [1000, 1000, 1000]

#E  ν

material_properties = [(29000,0.30, σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)]


# steel_law = MaterialModels.steel(σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)


# using Dierckx
# E_tan_spline = Spline1D(ϵ_s, E_tan)

# cdm = spl(ϵ_s)

#kx ky kϕ hx hy
springs = [0.0, 0.0, 0.0, 0.0, 0.0]

#member information
#L dL SectionProperties MaterialProperties Springs
member_definitions = [(13.0*12, 1.0,1,1,1)]

num_nodes = Int(member_definitions[1][1]/member_definitions[1][2]) + 1

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0, 13.0*12]

#imperfection
L = member_definitions[1][1]
dL = member_definitions[1][2]
z = 0: dL :L
Δo = L / 1000 * sin.(π*z/L) 

uo = Δo  #in the weak axis
vo = zeros(Float64, num_nodes)
ϕo = zeros(Float64, num_nodes)

imperfections = [uo, vo, ϕo]

#Calculate weak axis flexural buckling load
E = material_properties[1][1]
Iy = shape_info.Iy
Pcry = π^2*E*Iy/L^2   #weak axis buckling load

#calculate yield load
Py = 50*shape_info.A

#inelasticity
inelasticity_flag = 0

P = 0.88 * Py

loads = P * ones(num_nodes)


strain_limit = ϵy
# solution_controls = Column.SolutionControls("follow stress-strain curve", strain_limit)
# solution_controls = Column.SolutionControls("first yield", strain_limit)
solution_controls = Column.SolutionControls(" ", strain_limit)

#define column information
properties = Column.initialize(member_definitions, section_properties, section_geometry, material_properties, loads, springs, end_boundary_conditions, supports, imperfections, solution_controls)


#solve for column deformation
u, v, ϕ = Column.solve(properties, solution_controls)


# solution_controls = Column.SolutionControls(" ", strain_limit)
# u, v, ϕ = Column.solve(properties, solution_controls)


Myy = InternalForces.moment(properties.z, properties.dm, -u, properties.E, properties.Iy)



function calculate_strain(P, Myy, A, Iyy, cx, cxi, E)

    num_cells = length(cxi)
    ϵ_axial = zeros(Float64, num_cells)
    ϵ_flexure = zeros(Float64, num_cells)
    ϵ_total = zeros(Float64, num_cells)
 
    for i = 1:1:num_cells
          ϵ_axial[i] = -P/(A * E) 
          ϵ_flexure[i] = Myy/(E * Iyy) * (cxi[i] - cx)
          ϵ_total[i] = ϵ_axial[i] + ϵ_flexure[i]
 
    end
 
    return ϵ_axial, ϵ_flexure, ϵ_total
   
 end
 
 
 function update_cell_area(E_tan_spline, ϵ, Aio, E)

    num_cells = length(Aio)
    Ai = zeros(Float64, num_cells)
    E_tan = zeros(Float64, num_cells)
 
    for i = 1:num_cells
 
        # order = 1
        # x = (abs(ϵ[i]) - 2 * abs(ϵ[i])*1/1000):abs(ϵ[i])*1/1000:(abs(ϵ[i]) + 2 * abs(ϵ[i])*1/1000)
        # x0 = x[1]
        # sig = stress_strain_curve.(x)
        
        # stencil = Mesh.calculate_weights(order, x0, x)
 
        E_tan[i] = E_tan_spline(abs(ϵ[i]))
 
        Ai[i] = Aio[i] * E_tan[i]/E
 
    end
 
    return Ai, E_tan
 
 end

#    for i = 1:num_sections

    i = 78
      #calculate strain
      ϵ_axial, ϵ_flexure, ϵ_total = calculate_strain(properties.P[i], Myy[i], properties.A[i], properties.Iy[i], properties.xc[i], properties.cxi, properties.E[i])

      #update cross-section properties 
      properties.Ai, E_tan_track = update_cell_area(E_tan_spline, ϵ_total, properties.Aio, properties.E[i])
      properties.xc[i], properties.yc[i] = CrossSection.centroid_from_cells(properties.Ai, properties.cxi, properties.cyi)
      # properties.A[i] = CrossSection.area_from_cells(properties.Ai)
      properties.Iy[i] = CrossSection.moment_of_inertia_from_cells(properties.Ai, properties.cxi, properties.xc[i])

      #catch numerical issues I noticed, need to dig into this more
      if properties.A[i] <= 0.0
         # properties.A[i] = 0.0
         properties.Iy[i] = 0.0
      end

#    end



#if column yields, set I = 0 at the location of first yield.



#calculate moment
Myy = InternalForces.moment(properties.z, properties.dm, -u, properties.E, properties.Iy)

#baseline cell areas
Aio = deepcopy(Ai)

function calculate_strain(P, Myy, A, Iyy, cx, cxi, E)

    num_cells = length(cxi)
    ϵ_axial = zeros(Float64, num_cells)
    ϵ_flexure = zeros(Float64, num_cells)
    ϵ_total = zeros(Float64, num_cells)

    for i = 1:1:num_cells
        ϵ_axial[i] = -P/(A * E) 
        ϵ_flexure[i] = Myy/(E * Iyy) * (cxi[i] - cx)
        ϵ_total[i] = ϵ_axial[i] + ϵ_flexure[i]

    end

    return ϵ_axial, ϵ_flexure, ϵ_total

end


function update_cross_section_properties(Ai, cxi, cyi, P, Myy, E)

    cx = zeros(Float64, 10)
    Iyy = zeros(Float64, 10)

    #initialize 
    cx[1], cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)
    A = CrossSection.area_from_cells(Ai)
    Iyy[1] = CrossSection.moment_of_inertia_from_cells(Ai, cxi, cx[1])

    for i = 2:10
    
        #calculate strain
        ϵ_axial, ϵ_flexure, ϵ_total = calculate_strain(P, Myy, A, Iyy[i-1], cx[i-1], cxi, E)

        #update cross-section properties 
        Ai = update_cell_area(stress_strain_curve, ϵ_total, Aio, E)
        cx[i], cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)
        A = CrossSection.area_from_cells(Ai)
        Iyy[i] = CrossSection.moment_of_inertia_from_cells(Ai, cxi, cx[i])

        change_in_cx = abs((cx[i] - cx[i-1])/cx[i-1])

        if change_in_cx < 0.01

            return cx, Iyy

        else

            println("Number of maximum iterations exceeded.")
            return cx, Iyy

        end

    end

end


#calculate tangent elastic modulus
    



# cx, Iyy = update_cross_section_properties(Ai, cxi, cyi, P, Myy[78], Es)

 

# cx_it, Iyy_it = update_cross_section_properties(Ai, cxi, cyi, P, Myy[78], E)


cx = zeros(Float64, 10)
Iyy = zeros(Float64, 10)
A = zeros(Float64, 10)

#initialize 
cx[1], cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)
A[1] = CrossSection.area_from_cells(Ai)
Iyy[1] = CrossSection.moment_of_inertia_from_cells(Ai, cxi, cx[1])

for i = 2:10


    #calculate strain
    ϵ_axial, ϵ_flexure, ϵ_total = calculate_strain(P, Myy[78], A[i-1], Iyy[i-1], cx[i-1], cxi, Es)

    #update cross-section properties 
    Ai, E_tan = update_cell_area(stress_strain_curve, ϵ_total, Aio, Es)
    cx[i], cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)
    A[i] = CrossSection.area_from_cells(Ai)
    Iyy[i] = CrossSection.moment_of_inertia_from_cells(Ai, cxi, cx[i])


    ϵ_axial_i, ϵ_flexure_i, ϵ_total_i = calculate_strain(P, Myy[78], A[i], Iyy[i], cx[i], cxi, Es)

    ϵ_total_i_max = maximum(abs.(ϵ_total_i))   #new
    ϵ_total_max = maximum(abs.(ϵ_total))  #old

    if i > 2

        change_in_ϵ_total = abs((ϵ_total_i_max - ϵ_total_max)/ϵ_total_max)

        if change_in_ϵ_total < 0.01

            return cx, Iyy
    
        elseif i==10
    
            println("Number of maximum iterations exceeded.")
            return cx, Iyy
    
        end


    end

end
    # ϵ = 0.0001
    
    # order = 1
    # x = (abs(ϵ) - abs(ϵ)*1/1000):abs(ϵ)*1/1000:(abs(ϵ) + abs(ϵ)*1/1000)
    # x = [0.0, 0.00001, 0.00002]
    # x0 = x[1]
    # sig = stress_strain_curve.(x)
    # sig = 29000.0 .* x
    
    # stencil = Mesh.calculate_weights(order, x0, x)

    # E_tan = sum(stencil .* sig)


# end

Myy = InternalForces.moment(properties.z, properties.dm, -u, properties.E, properties.Iy)
ϵ_axial, ϵ_flexure, ϵ_total = calculate_strain(P, Myy[78], A[i-1], Iyy[i-1], cx[i-1], cxi, Es)

