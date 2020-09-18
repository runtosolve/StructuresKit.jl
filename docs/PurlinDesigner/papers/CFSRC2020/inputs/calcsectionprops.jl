using StructuresKit
using CSV, DataFrames


#generate centerline geometry of each specimen to calculate section properties
#use Python sectionproperties for this...


#read in cross-section dimensions from Gao and Moen (2013)
path = @__DIR__
filename = string(path, "/Gao_Moen_2013_Table1", ".csv")
dims = CSV.read(filename, DataFrame, header=false)

dims[1,13]

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

    testnum = string(i+1)


    filename = string("/Users/crismoen/RunToSolve/Papers/Moen_2020_CFSRC_PurlinDesigner", "/section_export", "/Test", testnum, "nodes", ".csv")
    CSV.write(filename,  DataFrame(node), header=false)

    filename = string("/Users/crismoen/RunToSolve/Papers/Moen_2020_CFSRC_PurlinDesigner", "/section_export", "/Test", testnum, "elements", ".csv")
    CSV.write(filename,  DataFrame(elem), header=false)

end




using Plots
plot(node[:,2], node[:, 3], markershape = :circle)





kϕ = 1339  #N-mm/rad/mm   c=34 mm closest to test

#need
#section properties
#local buckling
