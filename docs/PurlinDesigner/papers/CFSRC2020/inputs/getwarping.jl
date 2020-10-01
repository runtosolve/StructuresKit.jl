using StructuresKit
using CSV, DataFrames


#read in cross-section dimensions from Gao and Moen (2013)
path = @__DIR__
filename = string(path, "/Gao_Moen_2013_Table1", ".csv")
dims = CSV.read(filename, DataFrame, header=false)

Wn = zeros(Float64, (49, 2))

for i in 1:49


    #bring in cross-section dimensions
    Bc = dims[i,1]
    Dc = dims[i,2]
    θc = dims[i,4]

    Bt = dims[i,5]
    Dt = dims[i,6]
    θt = dims[i,8]

    H = dims[i,9]
    r = dims[i,10]
    t = dims[i,11]


    #calculate centerline dimensions

    CorZ = dims[i,13]+1
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


    prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,Bc,Bt,Dc,Dt,r,r,r,r,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)
    A, xc, yc, Ixc, Iyc, Ixyc, Imax, Imin, Th_p, Cw, J, Xs, Ys, w, Bx, By, B1, B2 = CrossSection.warp(node, elem)

    Wn[i, 1] = maximum(abs.(w))

end



filename = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/Wn.csv"
CSV.write(filename,  DataFrame(Wn), header=false)


############
