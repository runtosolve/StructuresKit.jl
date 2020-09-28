using CSV, DataFrames

testnum = 49

filename = string("/Users/crismoen/RunToSolve/Papers/Moen_2020_CFSRC_PurlinDesigner", "/section_export", "/Test", testnum, "nodes", ".csv")
nodes = CSV.read(filename, DataFrame, header=false)

path = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/inputs/"

filename = string(path, "julia_sectionprops.csv")
sectprops = CSV.read(filename, DataFrame, header=false)

filename = string(path, "/Gao_Moen_2013_Table1", ".csv")
dims = CSV.read(filename, DataFrame, header=false)


I = sectprops[testnum-1, 1]

Ma = 1.0
Mb = 2.0

Yc = sectprops[testnum-1,7]   #from bottom fiber

#calculate stress
σa = Ma.*(nodes[:,3] .- Yc)./I
σb = Mb.*(nodes[:,3] .- Yc)./I

dσ = σb .- σa

plot(nodes[:,2], nodes[:,3], dσ)
plot!(nodes[:,2], nodes[:,3], zeros(length(q)))


dx = 1.0


q = zeros(Float64, numnodes)
s = zeros(Float64, numnodes)

for i = 1:numelem

    x1 = nodes[i, 2]
    y1 = nodes[i, 3]
    x2 = nodes[i+1, 2]
    y2 = nodes[i+1, 3]

    dS = norm(x2-x1, y2-y1)
    dA = dS*t

    dq = (dσ[i]*dA + (1/2*(dσ[i]+dσ[i+1]) * dA))/dx

    q[i+1] = q[i] + dq
    s[i+1] = s[i] + dS

end


plot(nodes[:,2], nodes[:,3], q)
plot!(nodes[:,2], nodes[:,3], zeros(length(q)))

plot(s, q)

t = dims[testnum-1,11]



numnodes = length(nodes[:,1])
numelem = numnodes - 1

q = zeros(Float64, numnodes)
s = zeros(Float64, numnodes)

using LinearAlgebra
using Statistics

for i = 1:numelem

    x1 = nodes[i, 2]
    y1 = nodes[i, 3]
    x2 = nodes[i+1, 2]
    y2 = nodes[i+1, 3]

    dS = norm(x2-x1, y2-y1)
    y = Yc - mean([y1, y2])
    dq = V/I*y*t*dS

    q[i+1] = q[i] + dq
    s[i+1] = s[i] + dS

end

using Plots
plot(s, q)

plot(nodes[:,2], nodes[:,3], q)
plot!(nodes[:,2], nodes[:,3], zeros(length(q)))
