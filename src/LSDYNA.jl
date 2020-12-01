module LSDYNA

using Printf

export cross_section_extruder


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

    @printf(io, "ExecuteCommand(\"elgenerate shell curvedrag 1 1 %s %s 0 0 0 1 0 0 0",string(L), string(Int(L/dL)))

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

end #module