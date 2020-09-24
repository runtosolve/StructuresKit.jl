

testnum = 2

filename = string("/Users/crismoen/RunToSolve/Papers/Moen_2020_CFSRC_PurlinDesigner", "/section_export", "/Test", testnum, "nodes", ".csv")
CSV.read(filename,  DataFrame(node), header=false)

path = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/inputs/"

filename = string(path, "julia_sectionprops.csv")
sectprops = CSV.read(filename, DataFrame, header=false)


I = 1.0
V =
