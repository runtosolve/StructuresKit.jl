using StructuresKit



#A Ix Iy J Cw xo yo
SectionProperties = [(272.0,363370.0,64100.0, 188.0, 122720891.0, -32.59, 0.0)]


#E  ν
MaterialProperties = [(200,0.30)]


#kx ky kϕ hx hy
Springs = [(0.0, 0.0, 0.0, 0.0, 0.0)]

#member information
#L dL SectionProperties MaterialProperties Springs
MemberDefinitions = [(2438.0,2438.0/12,1,1,1)]

#P qx qy ax ay
Loads = [(1000*ones(13)),(2*ones(13)), (0.0*ones(13)),(0.0*ones(13)),(0.0*ones(13))]


#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
EndBoundaryConditions = [1 1]

#supports
#location where u=v=ϕ=0
Supports = [0.0 2438.0]


u, v, ϕ, properties = BeamColumn.solve(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports)
