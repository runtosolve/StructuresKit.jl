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

k = zeros(Float64, (14,2))


for i = 1:14

        ASDorLRFD=2    #ASDorLRFD=0 for ASD, =1 for LRFD

        MemberDefinitions = [(7468.0,7468.0/25,    1,1,1,1,1)]

        PurlinSpacing=(4140.0-2286.0)/2 + 2286.0/2;  #in.

        RoofSlope = 0.0;   #degrees

        #location where u=v=ϕ=0
        Supports = [0.0 7468.0]

        #end boundary conditions
        #type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
        EndBoundaryConditions = [1 1];

        #Each row below defines a set of section properties.  There are 4 rows because there are 4 cross-sections in this example -  8ZS2.25x070, 8ZS2.25x059, 8ZS2.25x070+ 8ZS2.25x059 at a splice, and 8ZS2.25x059+8ZS2.25x059 at a splice.

        #Ixy is negative here because the AISI D100 coordinate system has y pointing up however the PlautBeam.jl coordinate system has y pointing down.

        #Ix Iy Ixy J Cw Wn Mcrℓx Mcrℓy
        Ixx = sectprops[i,1]
        Iyy = sectprops[i,2]
        Ixy = sectprops[i,3]
        J = sectprops[i,4]
        Cw = sectprops[i,5]
        Wn = Wnall[i,1]
        Mcrlxx = Mcrlxxall[i,1]
        Mcrlyy = 99999999999

        SectionProperties = [(Ixx, Iyy, Ixy, J, Cw, Wn, Mcrlxx, Mcrlyy)];


        #t is base metal thickness
        #ho, b, d, and h are outside purlin depth, flange width, lip length, and web flat height
        #θ is lip angle from the horizontal
        #CorZ=0 for C, CorZ=1 for Z
        #This nomenclature is consistent with AISI S100-16.

        #t, ho, b, d, θ, CorZ, h
        bflange = dims[i,1]
        d = dims[i,2]
        θ = dims[i,4]
        ho = dims[i,9]
        r = dims[i,10]
        t = dims[i,11]
        h = ho - 2*(r+t)
        CorZ = dims[i, 13]


        # CrossSectionDimensions =
        # [(t, ho, b, d, θ , CorZ, h)];


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
        b = bflange - c

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

        #fastener spacing
        S = 305 #mm
        kp = 320 #N/mm   #Gao and Moen (2013) Fig. 9, 26 Gauge (t=0.46 mm)

        E=MaterialProperties[1][1]

        #calculate kϕ
        k[i, 1] = Connections.cfs_rot_screwfastened_k(b, c, S, t, kp, E, CorZ)


        t1 = 0.46   #26 gauge panel
        t2 = t
        E1 = E
        E2 = E
        Fu1 = Fy * 1.2
        Fu2 = Fy * 1.2
        Fss = 11100  #N   #12 screw
        D = 5.40  #  mm  #12 screw
        Ka, ψ, α, β, Ke = Connections.cfs_trans_screwfastened_k(t1, t2, E1, E2, Fss, Fu1, Fu2, D)

        k[i, 2] = Ke/S


end


filename = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/inputs/k.csv"
CSV.write(filename,  DataFrame(k), header=false)
