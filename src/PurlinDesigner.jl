module PurlinDesigner

#bring in local modules
using ..AISIS10016
using ..AISIS10024
using ..PlautBeam
using ..InternalForces
using ..BeamMesh

export lineStrength


#calculated section strengths along a purlin line
function sectionStrengths(ASDorLRFD, dz, dm, memberDefinitions, sectionProperties, crossSectionDimensions, materialProperties, loadLocation, bracingProperties, roofSlope)

    numberOfNodes = length(dz)+1

    #calculate strength limit state capacities

    #strong axis flexure, local-global interaction
    Fy = BeamMesh.propvector(memberDefinitions, dm, dz, materialProperties, 4, 3)
    Ixx = BeamMesh.propvector(memberDefinitions, dm, dz, sectionProperties, 3, 1)
    ho = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 7, 2)
    ycy = ho./2  #distance from neutral axis to outer fiber
    Sxx = Ixx./ycy
    Myxx = Fy.*Sxx
    Mnexx = Myxx
    Mcrℓxx = BeamMesh.propvector(memberDefinitions, dm, dz, sectionProperties, 3, 6)
    eMnℓxx =  AISIS10016.f321.(Mnexx, Mcrℓxx, ASDorLRFD)

    #weak axis flexure, local-global interaction
    Iyy = BeamMesh.propvector(memberDefinitions, dm, dz, sectionProperties, 3, 2)

    t = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 7, 1)
    b = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 3, 3)
    d = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 3, 4)
    θc = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 3, 5)

    #distance from neutral axis to outer fiber
    ycx = b.+d.*cos.(deg2rad.(θc)) .-t./2
    Syy = Iyy./ycx
    Myyy = Fy.*Syy
    Mneyy = Myyy
    Mcrℓyy = BeamMesh.propvector(memberDefinitions, dm, dz, sectionProperties, 3, 7)
    eMnℓyy = AISIS10016.f321.(Mneyy, Mcrℓyy, ASDorLRFD)

    #torsion
    Cw = BeamMesh.propvector(memberDefinitions, dm, dz, sectionProperties, 3, 5)
    Fy = BeamMesh.propvector(memberDefinitions, dm, dz, materialProperties, 4, 3)
    Wn = BeamMesh.propvector(memberDefinitions, dm, dz, sectionProperties, 3, 6)
    eBn = AISIS10024.h411.(Cw, Fy, Wn, ASDorLRFD)

    #distortional buckling
    CorZ = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 3, 6)
    CorZ = trunc.(Int, CorZ)
    E = BeamMesh.propvector(memberDefinitions, dm, dz, materialProperties, 4, 1)
    μ = BeamMesh.propvector(memberDefinitions, dm, dz, materialProperties, 4, 2)
    G= E./(2 .*(1 .+ μ))
    f1=1.0*ones(length(numberOfNodes))
    f2=-1.0*ones(length(numberOfNodes))
    Lm = BeamMesh.propvector(memberDefinitions, dm, dz, bracingProperties, 6, 3)

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

    eMnd = AISIS10016.f411.(My, Mcrd, ASDorLRFD)


    #Web shear strength
    a = BeamMesh.propvector(memberDefinitions, dm, dz, bracingProperties, 6, 4)
    h = BeamMesh.propvector(memberDefinitions, dm, dz, crossSectionDimensions, 7, 7)
    kv  = AISIS10016.g233.(a, h)
    Fcr = AISIS10016.g232.(E, μ, kv, h, t)
    Vcr = AISIS10016.g231.(h, t, Fcr)
    eVn = AISIS10016.g21.(E, h, t, Fy, Vcr, ASDorLRFD)


    #note nominal capacities are divided by ASD safety factor or
    #multiplied by LRFD resistance factor here
    return eMnℓxx, eMnℓyy, eBn, eMnd, eVn

end

#find flexure+torsion D/C
function bendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)

    #check bending + torsion interaction
    #ActionM1, ActionM2, ActionB, Interaction
    interactionCheck = AISIS10024.h42.(Mxx,Myy,B,eMnℓxx,eMnℓyy,eBn)

    actionMxx = [x[1] for x in interactionCheck]
    actionMyy = [x[2] for x in interactionCheck]
    actionB = [x[3] for x in interactionCheck]
    totalInteraction = [x[4] for x in interactionCheck]

    demandToCapacity = totalInteraction./1.15

    return actionMxx, actionMyy, actionB, totalInteraction, demandToCapacity

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
function lineDC(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

    z, u, v, ϕ, beamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad)

    Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
    Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

    BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = bendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
    distDemandToCapacity = distortionalDemandToCapacity(Mxx,eMnd)
    MVDemandToCapacity = bendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
    BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = biaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
    purlinDemandToCapacity=maximum([BTDemandToCapacity; distDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

    return purlinDemandToCapacity

end

#purlin line load that causes failure
function rootfinder(q, eps, residual, demandToCapacity, memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, loadAngle)

    #use bisection method of rootfinding to solve for purlin strength

    maxIterations=10  #hard code this for now

    for i=1:maxIterations

        newq=q/demandToCapacity
        q = q + (newq-q)/2

        qx = -q*sin(deg2rad(loadAngle))
        qy = q*cos(deg2rad(loadAngle))
        uniformLoad=(qx,qy)

        demandToCapacity=lineDC(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

        residual=1.0 - abs(demandToCapacity)

        if residual < eps
            return q
        end

    end

    error("Maximum rootfinding interation limit exceeded.  A failure load was not found.")

end

function lineStrength(ASDorLRFD, gravityOrUplift, memberDefinitions, sectionProperties, crossSectionDimensions, materialProperties, loadLocation, bracingProperties, roofSlope, endBoundaryConditions, supports)

    dz, z, dm = BeamMesh.define(memberDefinitions)

    #calculate expected capacities
    eMnℓxx, eMnℓyy, eBn, eMnd,  eVn  = sectionStrengths(ASDorLRFD, dz, dm, memberDefinitions, sectionProperties, crossSectionDimensions, materialProperties, loadLocation, bracingProperties, roofSlope)

    #use rootfinding to solve for expected strength

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
    demandToCapacity = lineDC(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

    #initialize residual for rootfinding
    residual=abs(1.0-demandToCapacity)
    eps=0.01  #residual tolerance, hard coded for now

    #this spits out the expected failure load
    q = rootfinder(q, eps, residual, demandToCapacity, memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, loadAngle)
    expectedPurlinLineStrength = q

    qx = -q*sin(deg2rad(loadAngle))
    qy = q*cos(deg2rad(loadAngle))
    uniformLoad = (qx,qy)

    #calculate all the final properties and deformations for the purlin line, at failure

    z, u, v, ϕ, beamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad)

    Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
    Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
    T = InternalForces.torsion(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.G, beamProperties.J, beamProperties.Cw)
    B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

    BTActionMxx, BTActionMyy, BTActionB, BTTotalInteraction, BTDemandToCapacity = bendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
    distDemandToCapacity = distortionalDemandToCapacity(Mxx,eMnd)
    MVDemandToCapacity = bendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
    BBActionP, BBActionMxx, BBActionMyy, BBTotalInteraction, BBDemandToCapacity = biaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
    demandToCapacity=maximum([BTDemandToCapacity; distDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

    return expectedPurlinLineStrength, z, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, Mxx, Myy, Vyy, T, B, BTActionMxx, BTActionMyy, BTActionB, BTTotalInteraction, BTDemandToCapacity, distDemandToCapacity, MVDemandToCapacity, BBActionP, BBActionMxx, BBActionMyy, BBTotalInteraction, BBDemandToCapacity, demandToCapacity

end

end #module
