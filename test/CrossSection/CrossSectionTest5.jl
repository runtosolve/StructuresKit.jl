using CSV, DataFrames

struct WShape

    d::Float64
    tw::Float64
    bf::Float64
    tf::Float64
    kdes::Float64
    k1::Float64
    n::Tuple

end


path = "/Users/crismoen/.julia/dev/StructuresKit/test/CrossSection/"
filename = string(path, "aisc-shapes-database-v15.0.csv")

shape_name = "W14X90"
n=(4, 4, 4, 4, 4)


data = CSV.read(filename, DataFrame, header=true)

shape_row = findfirst(==(shape_name), data.AISC_Manual_Label)



d = parse(Float64, data[shape_row, :d])
tw = parse(Float64, data[shape_row, :tw])
bf = parse(Float64, data[shape_row, :bf])
tf = parse(Float64, data[shape_row, :tf])
kdes = parse(Float64, data[shape_row, :kdes])

#get k1 from AISC table, it is in fraction format
k1 = data[shape_row, :k1]
index = findfirst("/", k1)

whole = Int(data[shape_row, :k1][1] - '0')

if isempty(index) == false
    top_fraction = parse(Float64, data[shape_row, :k1][index[1]-2:index[1]-1])
    bottom_fraction = parse(Float64, data[shape_row, :k1][index[1]+1:index[1]+2])
else
    top_fraction = 0.0
    bottom_fraction = 0.0
end

k1 = whole + top_fraction/bottom_fraction


section = WShape(d, tw, bf, tf, kdes,k1, n)

#from bottom of bottom flange, web centerline to left edge
xcoords = zeros(n[1]+1)
ycoords = zeros(n[1]+1)

flange_range = 0.0 : -section.bf / 2 / n[1] : -section.bf / 2
[xcoords[i] =  flange_range[i] for i in eachindex(flange_range)]
ycoords .= 0.0

#up along bottom flange thickness
flange_thickness_range = section.tf/n[2]:section.tf/n[2]:section.tf
xcoords = [xcoords; ones(n[2])*xcoords[end]]
ycoords = [ycoords; flange_thickness_range]

#over to fillet radius at bottom flange - web intersection

flange_flat = section.bf/2 - section.tw/2 - section.k1/2

inside_flange_range = (xcoords[end] + flange_flat/n[3]) : flange_flat/n[3] : (xcoords[end] + flange_flat)

xcoords = [xcoords; inside_flange_range]
ycoords = [ycoords; ones(n[3])*ycoords[end]]

#go around the fillet
radius = -xcoords[end] - section.tw
θ = (-π/2 + π/2/n[4]):π/2/n[4]: 0.0

xo = xcoords[end]
yo = ycoords[end] + radius

x_radius = xo .+ radius .* cos.(θ)
y_radius = yo .+ radius .* sin.(θ)

# plot(x_radius, y_radius, markershape = :o, linetype = :scatter)

xcoords = [xcoords; x_radius]
ycoords = [ycoords; y_radius]

#add web flat 
web_flat = section.d/2 - section.tf - radius

web_flat_range = (ycoords[end] + web_flat/n[5]): web_flat/n[5]: (ycoords[end] + web_flat)
xcoords = [xcoords; ones(n[5])*xcoords[end]]
ycoords = [ycoords; web_flat_range]

#mirror about horizontal axis
ycoords_horz_flip = ycoords .- ycoords[end]
ycoords_horz_flip = -ycoords_horz_flip
ycoords_horz_flip = ycoords_horz_flip .+ ycoords[end]

xcoords = [xcoords; reverse(xcoords)[2:end]]
ycoords = [ycoords; reverse(ycoords_horz_flip)[2:end]]

#mirror about vertical axis
xcoords_vert_flip = reverse(-xcoords)[2:end-1]

xcoords = [xcoords; xcoords_vert_flip]
ycoords = [ycoords; reverse(ycoords)[2:end-1]]


using Plots
plot(xcoords, ycoords, markershape = :o, linetype = :scatter, legend = false)


using StructuresKit

x = xcoords
y = ycoords
z = 0:1:120

coordinates = Mesh.extrude(xcoords, ycoords, z)

num_cross_section_nodes = length(x)
num_extruded_layers = length(z)
num_cross_section_elements = num_cross_section_nodes 


connectivity = Mesh.extruded_surface(num_cross_section_nodes, num_cross_section_elements, num_extruded_layers)

  #combine all the triangles









# coordinates, connectivity = Visualize.Mesh3D(cross_section_nodes, z, u, v, ϕ)

# #plot the mesh with Makie
using Makie
scene = poly(coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:false, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 4)



# plot(xcoords, ycoords, markershape = :o, linetype = :scatter, ylims=(0, 8))


#load CSV
#define cross-section geometry 
#extrude geometry
#define extruded facets
#define end facets 


