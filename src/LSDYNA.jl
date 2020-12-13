module LSDYNA

using Printf

export cross_section_extruder, cfile_cross_section_extruder


function cross_section_extruder(xcoords, ycoords, element_info, L, dL, path, scl_filename,lsdyna_filename)

    cd(path)

    io = open(scl_filename, "w")

    println(io,"define:")
    println(io,"void main(void)")
    println(io,"{")

    num_elem = size(element_info)[1]

    for i=1:num_elem
        @printf(io, "ExecuteCommand(\"line param  0 %f %f 0 %f %f\");\n", xcoords[i],ycoords[i],xcoords[i+1],ycoords[i+1])
    end

    println(io, "ExecuteCommand(\"elgenerate shelltype 6\");")
    println(io, "ExecuteCommand(\"genselect target occobject\");")
    println(io, "ExecuteCommand(\"occfilter add Edge\");")
    println(io, "ExecuteCommand(\"genselect whole\");")

    @printf(io, "ExecuteCommand(\"elgenerate shell curvedrag 1 1 %s %s 0 0 0 1 0 0 0",string(L), string(floor(L/dL)))

    for i=1:num_elem
        @printf(io," %s",string(i))
    end

    @printf(io, "\");\n")

    println(io, "ExecuteCommand(\"elgenerate accept\");")

    println(io, "ExecuteCommand(\"genselect whole\");")
    println(io, "ExecuteCommand(\"dupnode merge 0.001000\");")

    @printf(io, "ExecuteCommand(\"geomag del ")

    for i=1:num_elem
        @printf(io," %s",string(i, "e"))
    end

    @printf(io, "\");\n")

    println(io, string("ExecuteCommand(\"save keyword \\","\"",path,"/",lsdyna_filename,"\\","\"\");"))

    println(io,"}")
    println(io,"main();")

    close(io)

end


function cfile_cross_section_extruder(xcoords, ycoords, element_info, L, dL, cfile_filename,lsdyna_filename)

    # cd(path)

    io = open(cfile_filename, "w")

    num_elem = size(element_info)[1]

    for i=1:num_elem
        @printf(io, "line param  0 %f %f 0 %f %f\n", xcoords[i],ycoords[i],xcoords[i+1],ycoords[i+1])
    end

    println(io, "elgenerate shelltype 6")
    println(io, "genselect target occobject")
    println(io, "occfilter add Edge")
    println(io, "genselect whole")

    @printf(io, "elgenerate shell curvedrag 1 1 %s %s 0 0 0 1 0 0 0",string(L), string(floor(L/dL)))

    for i=1:num_elem
        @printf(io," %s",string(i))

        if i==num_elem
            @printf(io," %s\n",string(i))
        end
    end

    println(io, "elgenerate accept")

    # println(io, "genselect whole")
    # println(io, "dupnode merge 0.4000")

    @printf(io, "geomag del ")

    for i=1:num_elem
        @printf(io," %s",string(i, "e"))
        
        if i==num_elem
            @printf(io," %s\n",string(i, "e"))
        end
    end

  
    println(io, "dupnode open 1")
    @printf(io, "dupnode showdup %f\n", dL/100)
    println(io, "dupnode select 1")
    println(io, "genselect target node")
    println(io, "genselect node add part 1/0")
    @printf(io, "dupnode merge %f\n", dL/100)  #this dL/100 is a little tricky, might merge nodes not intended to be merged for a very fine mesh
    println(io, "genselect clear")
    println(io, "dupnode accept")

    # kfile_path="/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/Geometrica/cladding_study_2020/model_generation/all/"

    # lsdyna_filename=kfile_path*"h_purlin_"*string(i)*"_"*string(i)*".k"

    println(io, string("save keyword \"", lsdyna_filename, "\""))

    close(io)

end



end #module