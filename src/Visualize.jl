module Visualize

export Mesh3D


function Mesh3D(CrossSectionNodes, z, u, v, ϕ)

    #define planar cross-section geometry
    x=view(CrossSectionNodes,:,1)
    y=view(CrossSectionNodes,:,2)

    # write deformed cross-section in 3D coordinates, section by section along the z axis
    for i in 1:length(z)
        CrossSectionShape=hcat(x.+u[i]+x.*cos.(ϕ[i]).-y.*sin.(ϕ[i]),y.+v[i]+x.*sin.(ϕ[i]).+y.*cos.(ϕ[i]),ones(length(x))*z[i])
        if i==1
            global coordinates=CrossSectionShape   #start with first cross-section
        else
            coordinates=vcat(coordinates, CrossSectionShape)  #add following cross-sections along z
        end
    end

    # apply triangular mesh to 3D coordinates
    NumCrossSectionNodes=length(x)
    NumZNodes=length(z)

    #uppoly are triangles pointing towards +z, downpoly are triangles pointing towards -z
    uppoly=[[i i+1 i+1+NumCrossSectionNodes].+NumCrossSectionNodes*(j-1) for i=1:NumCrossSectionNodes-1 for j=1:NumZNodes-1]
    downpoly=[[i i+NumCrossSectionNodes i+NumCrossSectionNodes+1].+NumCrossSectionNodes*(j-1) for i=1:NumCrossSectionNodes-1 for j=1:NumZNodes-1]

    connectivity=vcat(uppoly,downpoly)  #combine all the triangles
    connectivity=permutedims(reshape(hcat(connectivity...), (length(connectivity[1]), length(connectivity))))  #convert into format to be read by Makie

    return coordinates, connectivity   #these are ready for Makie scenes, e.g., poly()

end

end #module
