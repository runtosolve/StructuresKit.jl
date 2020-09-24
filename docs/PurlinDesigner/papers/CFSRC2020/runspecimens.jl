using StructuresKit
using CSV, DataFrames

#Import Gao and Moen (2013) specimen section properties and Mcrl.

path = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/inputs/"

filename = string(path, "julia_sectionprops.csv")
sectprops = CSV.read(filename, DataFrame, header=false)

filename = string(path, "julia_Mcrl_xx.csv")
Mcrlxxall = CSV.read(filename, DataFrame, header=false)

filename = string(path, "Wn.csv")
Wnall = CSV.read(filename, DataFrame, header=false)

filename = string(path, "/Gao_Moen_2013_Table1", ".csv")
dims = CSV.read(filename, DataFrame, header=false)


using StructuresKit

FailurePressure = zeros(Float64, 49)

index = [5; 6; 7; 8; 11; 12]


# for j = 1:length(index)


        # i = index[j]

        i  = 49

        ASDorLRFD=2    #ASDorLRFD=0 for ASD, =1 for LRFD

        MemberDefinitions = [(7468.0,7468.0/100,    1,1,1,1,1)]

        PurlinSpacing=(4140.0-2286.0)/2 + 2286.0/2;  #in.

        RoofSlope = 0.0;   #degrees

        #location where u=v=ϕ=0
        Supports = [0.0 7468.0]

        #end boundary conditions
        #type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
        EndBoundaryConditions = [1 1];

        #Each row below defines a set of section properties.  There are 4 rows because there are 4 cross-sections in this example -  8ZS2.25x070, 8ZS2.25x059, 8ZS2.25x070+ 8ZS2.25x059 at a splice, and 8ZS2.25x059+8ZS2.25x059 at a splice.

        #Ixy is negative here because the AISI D100 coordinate system has y pointing up however the PlautBeam.jl coordinate system has y pointing down.

        #Ix Iy Ixy J Cw Wn Mcrℓx Mcrℓy Bcrℓ
        Ixx = sectprops[i,1]
        Iyy = sectprops[i,2]
        Ixy = -sectprops[i,3]
        J = sectprops[i,4]
        Cw = sectprops[i,5]
        Wn = Wnall[i,1]
        Mcrlxx = Mcrlxxall[i,1]
        Mcrlyy = 99999999999.0
        Bcrl = 99999999999.0

        SectionProperties = [(Ixx, Iyy, Ixy, J, Cw, Wn, Mcrlxx, Mcrlyy, Bcrl)];


        #t is base metal thickness
        #ho, b, d, and h are outside purlin depth, flange width, lip length, and web flat height
        #θ is lip angle from the horizontal
        #CorZ=0 for C, CorZ=1 for Z
        #This nomenclature is consistent with AISI S100-16.

        #t, ho, b, d, θ, CorZ, h
        b = dims[i,1]
        d = dims[i,2]
        θ = dims[i,4]
        ho = dims[i,9]
        r = dims[i,10]
        t = dims[i,11]
        h = ho - 2*(r+t)
        CorZ = dims[i, 13]


        CrossSectionDimensions =
        [(t, ho, b, d, θ , CorZ, h)];


        #ax ay


        #location of centroid
        Xc = sectprops[i,6]   #from centerline web
        Yc = sectprops[i,7]   #from bottom fiber

        #distance from centroid to shear center
        Xs = sectprops[i,8]
        Ys = sectprops[i,9]

        #location of shear center
        Xs_loc = Xc+Xs
        Ys_loc = Yc+Ys

        c = dims[i, 14]

        #for a C
        if CorZ == 0
                ax = abs(Xs_loc) - t/2 + b - c
        elseif CorZ == 1  #for a Z
                ax = abs(Xs_loc) - t/2 + c
        end

        ay = ho - abs(Ys_loc)  #top of top flange

        LoadLocation = [(ax, ay)];


        #E  ν  Fy
        Fy = dims[i, 12]
        MaterialProperties = [(203395,0.30, Fy)];



        #Bracing stiffness
        #Lm is the spacing between bracing that restrains distortional buckling
        #a is the web shear stiffener spacing, assumed equal to the span length here since none are provided
        #kx  kϕ  Lm  a

        kϕ = dims[i,15]*1.0
        kx = dims[i,16]*1.0


        BracingProperties = [(kx,kϕ, 7468.0, 7468.0)];



        GravityOrUplift=1   #GravityOrUplift=0 for gravity loading






        #calculate shear flow for Cee
        #find deformation
        #find demand moment
        #find capacity
        #add to interaction equations




#         eqn, z, deformation, strengths, forces, interactions, demand_to_capacity = PurlinDesigner.lineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)
#
#         FailurePressure[i]=eqn/PurlinSpacing*10^6   #N/m^2
#
#
#         UniformLoad=[0.0, eqn]
#
#         z, u, v, ϕ, beamProperties = PlautBeam.solve(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad);
#
#
#


#for uplift loading, consider the free flange
#calculate section properties

    #bring in cross-section dimensions
    Bc = dims[i,1]
    Dc = dims[i,2]
    θc = dims[i,4]

    Bt = dims[i,5]
    Dt = dims[i,6]
    θt = dims[i,8]

    H = dims[i,9]
    r = 0.0
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

    if CorZ == 2
        node[:, 2] = -node[:, 2]
    end

    if CorZ == 1
        index_xo = findall(x-> x==0.0, node[:, 2])
        index_yo = findall(y-> y==0.0, node[:, 3])
        index_o = intersect(index_x, index_y)
        index = 1:index_o[1]

elseif CorZ == 2
        index = findall(x->x<0, node[:, 2])
    end


    nodeflange = node[index,:]

      plot(nodeflange[:,2], nodeflange[:,3])

    elemflange = elem[index[1:end-1],:]


    A, xc, yc, Ixc, Iyc, Ixyc, Imax, Imin, Th_p, Cw, J, Xs, Ys, w, Bx, By, B1, B2 = CrossSection.CUFSMproperties(nodeflange, elemflange)

    Icantilever = 1/12*t^3   #length^4/length for distributed spring

    kx = 3*E*Icantilever/H^3
    kϕ = E*Icantilever/H




    #A Ix Iy J Cw xo yo
    SectionProperties = [(A,Ix,64100.0, 188.0, 122720891.0, -32.59, 0.0)]


    #kx ky kϕ hx hy
    Springs = [(0.0, 0.0, 0.0, 0.0, 0.0)]

    #member information
    #L dL SectionProperties MaterialProperties Springs
    MemberDefinitions = [(2438.0,2438.0/12,1,1,1)]

    #P qx qy ax ay
    Loads = [(10000*ones(13)),(2*ones(13)), (0.0*ones(13)),(0.0*ones(13)),(-1.0*ones(13))]


    #end boundary conditions
    #type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
    EndBoundaryConditions = [1 1]

    #supports
    #location where u=v=ϕ=0
    Supports = [0.0 2438.0]


   BeamColumn.solve(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports)



# # end
#
#
# using Plots
# plot(z, u)
#
# plot(z, forces.Myy)
#
#
# using DelimitedFiles
# filename = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/failurepressure.csv"
# writedlm(filename, FailurePressure)
#
# hall = dims[:,9]
# tall = dims[:,11]
#
# Ptest = dims[:, 17]
#
# using Plots
# index = [5; 6; 7; 8; 11; ;12]
# plot(index, -Ptest[index]./FailurePressure[index], seriestype = :scatter)
#
#
#
#
# plot(z, deformation.ϕ)
# plot(z, forces.T)
# plot!(z, strengths.eMnℓyy)
#
# plot(z, interactions.BTB)
#
# plot(z, demand_to_capacity.BT)
# plot(z, demand_to_capacity.dist)
#
# plot(z, )
#
# :BTMxx, :BTMyy, :BTB, :BTTotal, :BBP, :BBMxx, :BBMyy, :BBTotal
#
# :BT, :dist, :MV, :BB, :envelope
#
#
#
#
#
#
#
# println("ASD expected gravity roof capacity = ",round(FailurePressure,digits=1), " psf")
