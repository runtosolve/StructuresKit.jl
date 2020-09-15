using JuliaFEM
using JuliaFEM.Preprocess
add_elements! = JuliaFEM.add_elements!

# https://mechanicalc.com/reference/cross-sections

mesh = abaqus_read_mesh(joinpath("cross_section_properties", "beam.inp"))

volume_elements = create_elements(mesh, "BEAM")
left_elements = create_surface_elements(mesh, "LEFT")

info("Number of volume elements: ", length(volume_elements))
info("Number of surface elements in LEFT", length(left_elements))

update!(volume_elements, "density", 5.0)
volume = 0.0
mass = 0.0

# Calculate volume and mass of body

time = 0.0
for element in volume_elements
    for ip in get_integration_points(element)
        detJ = element(ip, time, Val{:detJ})
        volume += ip.weight * detJ
        density = element("density", ip, time)
        mass += ip.weight * density * detJ
    end
end

info("Volume of mesh: ", round(volume, 2))
info("Mass of mesh: ", round(mass, 2))

b = 10
h = 20
L = 100

using Base.Test
@test isapprox(volume, b*h*L)
@test isapprox(mass, 5*volume)

# Calculate first moment of inertia with respect to origin

Qy = 0.0
Qz = 0.0
area = 0.0
for element in left_elements
    for ip in get_integration_points(element)
        detJ = element(ip, time, Val{:detJ})
        x, y, z = element("geometry", ip, time)
        Qy += ip.weight * z * detJ
        Qz += ip.weight * y * detJ
        area += ip.weight * detJ
    end
end

info("Qy: ", round(Qy, 2))
info("Qz: ", round(Qz, 2))

# Calculate centroid

yc = Qz/area
zc = Qy/area

info("yc: ", round(yc, 2))
info("zc: ", round(zc, 2))

@test isapprox(yc, b/2)
@test isapprox(zc, h/2)

# Calculate second moment of inertia with respect to centroid

Iy = 0.0
Iz = 0.0
for element in left_elements
    for ip in get_integration_points(element, 1)
        detJ = element(ip, time, Val{:detJ})
        x, y, z = element("geometry", ip, time)
        Iy += ip.weight * (y-yc)^2 * detJ
        Iz += ip.weight * (z-zc)^2 * detJ
        area += ip.weight * detJ
    end
end

info("Iy: ", round(Iy, 2))
info("Iz: ", round(Iz, 2))

@test isapprox(Iy, h*b^3/12)
@test isapprox(Iz, b*h^3/12)
