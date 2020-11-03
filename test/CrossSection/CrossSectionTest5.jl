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
θ = π/2/n[4]:π/2/n[4]:π/2
x_radius = radius .* sin.(θ)
y_radius = radius .* cos.(θ)

#load CSV
#define cross-section geometry 
#extrude geometry
#define extruded facets
#define end facets 



