module PurlinDesigner

#bring in local modules
using ..AISIS10016
using ..AISIS10024
using ..Beam
using ..InternalForces
using ..Mesh
using ..BeamColumn
using ..CrossSection

export line_strength, free_flange_define


#calculated section strengths along a purlin line
function sectionStrengths(ASDorLRFD, dz, dm, memberDefinitions, sectionProperties, crossSectionDimensions, materialProperties, loadLocation, bracingProperties, roofSlope, FlangeProperties)

    numberOfNodes = length(dz)+1

    #calculate strength limit state capacities

    #strong axis flexure, local-global interaction
    Fy = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, materialProperties, 4, 3)
    Ixx = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 1)
    ho = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 7, 2)
    ycy = ho./2  #distance from neutral axis to outer fiber
    Sxx = Ixx./ycy
    Myxx = Fy.*Sxx
    Mnexx = Myxx
    Mcrℓxx = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 7)

    Mnℓxx = zeros(Float64, numberOfNodes)
    eMnℓxx = zeros(Float64, numberOfNodes)
    for i in eachindex(Mcrℓxx)
        Mnℓxx[i], eMnℓxx[i] =  AISIS10016.f321(Mnexx[i], Mcrℓxx[i], ASDorLRFD)
    end

    #weak axis flexure, local-global interaction
    Iyy = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 2)

    t = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 7, 1)
    b = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 3, 3)
    d = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 3, 4)
    θc = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 3, 5)

    #distance from neutral axis to outer fiber
    ycx = b.+d.*cos.(deg2rad.(θc)) .-t./2
    Syy = Iyy./ycx
    Myyy = Fy.*Syy
    Mneyy = Myyy
    Mcrℓyy = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 8)

    Mnℓyy = zeros(Float64, numberOfNodes)
    eMnℓyy = zeros(Float64, numberOfNodes)
    for i in eachindex(Mcrℓyy)
        Mnℓyy[i], eMnℓyy[i] = AISIS10016.f321(Mneyy[i], Mcrℓyy[i], ASDorLRFD)
    end

    #torsion
    Cw = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 5)
    Fy = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, materialProperties, 4, 3)
    Wn = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 6)
    Bcrℓ = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, sectionProperties, 3, 9)

    Bn = zeros(Float64, numberOfNodes)
    eBn = zeros(Float64, numberOfNodes)
    for i in eachindex(Bcrℓ)
        Bn[i], eBn[i] = AISIS10024.h411(Cw[i], Fy[i], Wn[i], Bcrℓ[i], ASDorLRFD)
    end

    #distortional buckling
    CorZ = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 3, 6)
    CorZ = trunc.(Int, CorZ)
    E = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, materialProperties, 4, 1)
    μ = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, materialProperties, 4, 2)
    G= E./(2 .*(1 .+ μ))
    f1=1.0*ones(length(numberOfNodes))
    f2=-1.0*ones(length(numberOfNodes))
    Lm = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, bracingProperties, 6, 3)

    #define moment gradient factor
    #assume it is 1
    M1=1.0*ones(length(numberOfNodes))
    M2=1.0*ones(length(numberOfNodes))
    curvatureSign = -1.0*ones(length(numberOfNodes))

    #assume clip stiffness does not help distortional buckling for a standing seam roof
    #clips are spaced too far apart
    kϕ=0.0*zeros(length(numberOfNodes))

    Sf=Sxx

    Mcrd=AISIS10016.app23331.(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, curvatureSign, Lm, kϕ, Sf)

    My=Fy.*Sf

    Mnd = zeros(Float64, numberOfNodes)
    eMnd = zeros(Float64, numberOfNodes)
    for i in eachindex(Mcrd)
        Mnd[i], eMnd[i] = AISIS10016.f411(My[i], Mcrd[i], ASDorLRFD)
    end


    #Web shear strength
    a = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, bracingProperties, 6, 4)
    h = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, crossSectionDimensions, 7, 7)
    kv  = AISIS10016.g233.(a, h)
    Fcr = AISIS10016.g232.(E, μ, kv, h, t)
    Vcr = AISIS10016.g231.(h, t, Fcr)

    Vn = zeros(Float64, numberOfNodes)
    eVn = zeros(Float64, numberOfNodes)
    for i in eachindex(Vcr)
        Vn[i], eVn[i] = AISIS10016.g21(E[i], h[i], t[i], Fy[i], Vcr[i], ASDorLRFD)
    end


    #define free flange properties
    Iyyf = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, FlangeProperties, 3, 2)
    xcf = Mesh.create_line_element_property_array(memberDefinitions, dm, dz, FlangeProperties, 3, 6)


    #Free flange strength
    #Assume yield criterion for now, no local buckling...
    Sycf_tension = Iyyf ./abs.(xcf)
    Sycf_comp = Iyyf ./((b .+d .*cos.(deg2rad.(θc)) .- t/2) .-abs.(xcf))
    Sycf = min.([Sycf_tension Sycf_comp])
    My_yyf = Fy.*Sycf
    Mcrℓ_yyf = 999999999999999.0  #no local buckling right now


    Mnℓyy_freeflange = zeros(Float64, numberOfNodes)
    eMnℓyy_freeflange = zeros(Float64, numberOfNodes)
    for i in eachindex(Mnℓyy_freeflange)
        Mnℓyy_freeflange[i], eMnℓyy_freeflange[i] = AISIS10016.f321(My_yyf[i], Mcrℓ_yyf, ASDorLRFD)
    end


    #note nominal capacities are divided by ASD safety factor or
    #multiplied by LRFD resistance factor here
    return eMnℓxx, eMnℓyy, eBn, eMnd, eVn, eMnℓyy_freeflange

end

#find flexure+torsion D/C
function bendingTorsionDemandToCapacity(Mxx, Myy, B, Myy_freeflange, eMnℓxx, eMnℓyy, eBn, eMnℓyy_freeflange)

    #check bending + torsion interaction
    #ActionM1, ActionM2, ActionB, Interaction
    interactionCheck = AISIS10024.h42.(Mxx, Myy, B, Myy_freeflange, eMnℓxx, eMnℓyy, eBn, eMnℓyy_freeflange)

    actionMxx = [x[1] for x in interactionCheck]
    actionMyy = [x[2] for x in interactionCheck]
    actionB = [x[3] for x in interactionCheck]
    actionMyy_freeflange = [x[4] for x in interactionCheck]
    totalInteraction = [x[5] for x in interactionCheck]

    demandToCapacity = totalInteraction./1.15

    return actionMxx, actionMyy, actionB, actionMyy_freeflange, totalInteraction, demandToCapacity

end

#find distortional buckling D/C
function distortionalDemandToCapacity(Mxx,eMnd)

    #check distortional buckling
    demandToCapacity=abs.(Mxx./eMnd)

end

#find flexure+shear D/C
function bendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)

    interaction = AISIS10016.h21.(Mxx, Vyy, eMnℓxx, eVn)

    demandToCapacity = interaction

    return demandToCapacity

end

#find biaxial bending D/C
function biaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)

    #no axial force for now
    Pbar=zeros(length(Mxx))
    Pa=ones(length(Mxx))

    interactionCheck = AISIS10016.h121.(Pbar, Mxx, Myy, Pa, eMnℓxx, eMnℓyy)

    actionP = [x[1] for x in interactionCheck]
    actionMxx = [x[2] for x in interactionCheck]
    actionMyy = [x[3] for x in interactionCheck]
    totalInteraction = [x[4] for x in interactionCheck]

    demandToCapacity = totalInteraction

    return actionP, actionMxx, actionMyy, totalInteraction, demandToCapacity

end

#purlin line demand to capacity
function lineDC(memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, eMnℓyy_freeflange, bridging)

    #Calculate purlin line deformation and demands.
    z, u, v, ϕ, beamProperties = Beam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad)
    Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
    Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

    #Calculate free flange behavior with shear flow and flexural-torsional buckling.
    FlangeProperties, Springs, Loads = PurlinDesigner.free_flange_define(memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, bracingProperties, uniformLoad[2], Mxx)
    uf, vf, ϕf, FreeFlangeProperties = BeamColumn.solve(memberDefinitions, FlangeProperties, materialProperties, Loads, Springs, endBoundaryConditions, bridging)
    Myyf = InternalForces.moment(z, beamProperties.dm, -uf, FreeFlangeProperties.E, FreeFlangeProperties.Iy)

    #Get limit state interactions and demand-to-capacities along the purlin line.
    BTActionMxx, BTActionMyy, BTActionB, BTActionMyy_freeflange, BTTotalInteraction, BTDemandToCapacity = bendingTorsionDemandToCapacity(Mxx, Myy, B, Myyf, eMnℓxx, eMnℓyy, eBn, eMnℓyy_freeflange)
    distDemandToCapacity = distortionalDemandToCapacity(Mxx,eMnd)
    MVDemandToCapacity = bendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
    BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = biaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)

    #Find worst case demand-to-capacity ratio along the purlin line.
    purlinDemandToCapacity=maximum([BTDemandToCapacity; distDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

    return purlinDemandToCapacity

end

#purlin line load that causes failure
function rootfinder(q, eps, residual, demandToCapacity, memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, loadAngle, eMnℓyy_freeflange, bridging)

    #use bisection method of rootfinding to solve for purlin strength

    maxIterations=10  #hard code this for now

    for i=1:maxIterations

        newq=q/demandToCapacity
        q = q + (newq-q)/2

        qx = -q*sin(deg2rad(loadAngle))
        qy = q*cos(deg2rad(loadAngle))
        uniformLoad=(qx,qy)


        demandToCapacity=lineDC(memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, eMnℓyy_freeflange, bridging)

        residual = 1.0 - abs(demandToCapacity)

        if abs(residual) < eps
            return q
        end

    end

    error("Maximum rootfinding interation limit exceeded.  A failure load was not found.")

end

function line_strength(ASDorLRFD, gravityOrUplift, memberDefinitions, sectionProperties, crossSectionDimensions, materialProperties, loadLocation, bracingProperties, roofSlope, endBoundaryConditions, supports, bridging)

    dz, z, dm = Mesh.define_line_element(memberDefinitions)
    numnodes = length(z)

    #Calculate free flange available strengths.
    FlangeProperties, Springs, Loads = PurlinDesigner.free_flange_define_line_element(memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, bracingProperties, 0.0, zeros(Float64, numnodes))


    #Calculate purlin available strengths.
    eMnℓxx, eMnℓyy, eBn, eMnd,  eVn, eMnℓyy_freeflange  = sectionStrengths(ASDorLRFD, dz, dm, memberDefinitions, sectionProperties, crossSectionDimensions, materialProperties, loadLocation, bracingProperties, roofSlope, FlangeProperties)

    #Use rootfinding to solve for purlin line available strength.

    #pick a small load
    if gravityOrUplift==0
      q=0.0001
      loadAngle=roofSlope
    elseif gravityOrUplift==1
      q=-0.0001
      loadAngle=0.0   #for uplift load is always perpendicular to purlin flange
    end

    #find DemandToCapacity for this small load
    qx = -q*sin(deg2rad(loadAngle))
    qy = q*cos(deg2rad(loadAngle))
    uniformLoad = (qx,qy)
    demandToCapacity = lineDC(memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, eMnℓyy_freeflange, bridging)

    #initialize residual for rootfinding
    residual=abs(1.0-demandToCapacity)
    eps=0.01  #residual tolerance, hard coded for now

    #this spits out the expected failure load
    eqn = rootfinder(q, eps, residual, demandToCapacity, memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, loadAngle, eMnℓyy_freeflange, bridging)

    qx = -eqn*sin(deg2rad(loadAngle))
    qy = eqn*cos(deg2rad(loadAngle))
    uniformLoad = (qx,qy)

    #calculate all the final properties and deformations for the purlin line, at failure

    z, u, v, ϕ, beamProperties = Beam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad)

    Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
    Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    T = InternalForces.torsion(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.G, beamProperties.J, beamProperties.Cw)
    B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

    #Calculate free flange behavior with shear flow and flexural-torsional buckling.
    FlangeProperties, Springs, Loads = PurlinDesigner.free_flange_define(memberDefinitions, sectionProperties, materialProperties, crossSectionDimensions, bracingProperties, uniformLoad[2], Mxx)
    uf, vf, ϕf, FreeFlangeProperties = BeamColumn.solve(memberDefinitions, FlangeProperties, materialProperties, Loads, Springs, endBoundaryConditions, bridging)

    Myyf = InternalForces.moment(z, beamProperties.dm, -uf, FreeFlangeProperties.E, FreeFlangeProperties.Iy)

    #Get limit state interactions and demand-to-capacities along the purlin line.
    BTActionMxx, BTActionMyy, BTActionB, BTActionMyy_freeflange, BTTotalInteraction, BTDemandToCapacity = bendingTorsionDemandToCapacity(Mxx, Myy, B, Myyf, eMnℓxx, eMnℓyy, eBn, eMnℓyy_freeflange)
    distDemandToCapacity = distortionalDemandToCapacity(Mxx,eMnd)
    MVDemandToCapacity = bendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
    BBActionP, BBActionMxx, BBActionMyy, BBTotalInteraction, BBDemandToCapacity = biaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)

    deformation = NamedTuple{(:u, :v, :ϕ, :uf, :vf, :ϕf)}((u, v, ϕ, uf, vf, ϕf))
    strengths = NamedTuple{(:eMnℓxx, :eMnℓyy, :eBn, :eMnd, :eVn, :eMnℓyy_freeflange)}((eMnℓxx, eMnℓyy, eBn, eMnd, eVn, eMnℓyy_freeflange))
    forces = NamedTuple{(:Mxx, :Myy, :Vyy, :T, :B, :Myyf)}((Mxx, Myy, Vyy, T, B, Myyf))
    interactions = NamedTuple{(:BTMxx, :BTMyy, :BTB, :BTMyy_freeflange, :BTTotal, :BBP, :BBMxx, :BBMyy, :BBTotal)}((BTActionMxx, BTActionMyy, BTActionB, BTActionMyy_freeflange, BTTotalInteraction, BBActionP, BBActionMxx, BBActionMyy, BBTotalInteraction))
    demand_to_capacity = NamedTuple{(:BT, :dist, :MV, :BB, :envelope)}((BTDemandToCapacity, distDemandToCapacity, MVDemandToCapacity, BBDemandToCapacity, demandToCapacity))


    return eqn, z, beamProperties, deformation, strengths, forces, interactions, demand_to_capacity

end


function free_flange_stiffness(t, E, H, kϕc)

    Icantilever = 1/12*t^3   #length^4/length for distributed spring


    #Use Eq. 16 from Gao and Moen (2013) https://ascelibrary.org/doi/abs/10.1061/(ASCE)ST.1943-541X.0000860
    #kxf = 3*E*Icantilever/H^3
    kxf = 1/(H^2/kϕc + (H^3/(3*E*Icantilever)))

    kϕf = E*Icantilever/H

    return kxf, kϕf

end


function free_flange_define(MemberDefinitions, SectionProperties, MaterialProperties, CrossSectionDimensions, BracingProperties, q, Mxx)


    dz, z, dm = Mesh.define_line_element(MemberDefinitions)

    numnodes = length(z)

    CorZ = CrossSectionDimensions[1][6] + 1
    H = CrossSectionDimensions[1][2]
    Bc = CrossSectionDimensions[1][3]
    Dc = CrossSectionDimensions[1][4]
    θc = CrossSectionDimensions[1][5]
    t = CrossSectionDimensions[1][1]

    r = 0.0
    kipin = 0
    center = 0
    n = 5

    node, elem = CrossSection.CZflange_template(CorZ,H,Bc,Bc,Dc,Dc,r,r,r,r,θc,θc,t,n,n,n,n,n,n,n,n,n,kipin,center)


    coords = node[:, 2:3]
    ends = elem[:, 2:4]

    Af,xcf,ycf,Ixf,Iyf,Ixyf,thetaf,I1f,I2f,Jf,xsf,ysf,Cwf,B1f,B2f,wnf = CrossSection.CUFSMsection_properties(coords, ends)



    FlangeProperties = [(Af, Ixf, Iyf, Jf, Cwf, xcf, ycf, xsf, ysf)]

    E = MaterialProperties[1][1]   #consider generalizing this someday

    kϕc = BracingProperties[1][2]

    kxf, kϕf = free_flange_stiffness(t, E, H, kϕc)

    #kx ky kϕ hx hy
    Springs = [(kxf*ones(numnodes)),(0.0*ones(numnodes)), (kϕf*ones(numnodes)),(0.0*ones(numnodes)),(0.0*ones(numnodes))]

    #approximate axial force in flange
    P = -Mxx ./ ((H -t) - 2 * abs(ycf))

    #There is shear flow in the free flange of a restrained C or Z.
    #Use Eq. 13b and 15 from Gao and Moen (2013) https://ascelibrary.org/doi/abs/10.1061/(ASCE)ST.1943-541X.0000860

    Ix = SectionProperties[1][1]
    c = CrossSectionDimensions[1][8]
    b = Bc - c

    if CorZ == 1  #C

        kH = ((Bc^2*t*H^2)/(4*Ix) + b)/H

        if xcf > 0   #if flange is oriented left to right from web
            kH = -kH
        end


    elseif CorZ == 2  #Z

        kH = (H*t*(Bc^2 + 2*Dc*Bc - (2*Dc^2*Bc)/H))/(4*Ix)

        if xcf > 0   #if flange is oriented left to right from web
            kH = -kH
        end

    end

    qx = q .* kH

    #P qx qy ax ay
    Loads = [P, (qx*ones(numnodes)), (0.0*ones(numnodes)),(0.0*ones(numnodes)),(ycf*ones(numnodes))]

    return FlangeProperties, Springs, Loads

end




end #module
