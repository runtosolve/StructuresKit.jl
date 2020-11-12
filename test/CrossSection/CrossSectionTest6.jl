using CSV, DataFrames


struct WShape

    d::Float64
    tw::Float64
    bf::Float64
    tf::Float64
    kdes::Float64
    k1::Float64
    A::Float64
    Ix::Float64
    Iy::Float64
    J::Float64
    Cw::Float64
    Zx::Float64
    Zy::Float64
    Wno::Float64


end

function AISC(shape_name)

    filename = string(@__DIR__, "/aisc-shapes-database-v15.0.csv")

    data = CSV.read(filename, DataFrame, header=true)

    shape_row = findfirst(==(shape_name), data.AISC_Manual_Label)

    section_type = data[shape_row, :Type]

    if section_type == "W"

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

        A = data[shape_row, :A]
        Ix = data[shape_row, :Ix]
        Iy = data[shape_row, :Iy]
        J = parse(Float64, data[shape_row, :J])
        Cw = parse(Float64, data[shape_row, :Cw])
        Zx = data[shape_row, :Zx]
        Zy = data[shape_row, :Zy]
        Wno = parse(Float64, data[shape_row, :Wno])

        section = WShape(d, tw, bf, tf, kdes,k1, A, Ix, Iy, J, Cw, Zx, Zy, Wno)

    end

    return section

end    

shape_name = "W14X90"


section = AISC(shape_name)