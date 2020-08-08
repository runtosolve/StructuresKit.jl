using StructuresKit

#362S162-54, Mcrd, no springs, no moment gradient, Lm>>Lcrd

CorZ=0
t=0.0566

#out to out dimensions
ho=3.625
b=1.625
d=0.5
θc=90.0

E=29500.0
μ=0.30
G=E/(2*(1+μ))
f1=1
f2=-1
M1=1.0
M2=1.0
Lm=100.0
kϕ=0.0

Sf=0.9085/(3.625/2)

CurvatureSign = -1.0

Mcrd = AISIS10016.app23331(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, CurvatureSign, Lm, kϕ, Sf)


#compare to CFSEI TechNote G100-08 (or G100-09?), Table 3
#https://cfsei.memberclicks.net/technical-notes
@test Mcrd ≈ 89.2*Sf atol=0.01*Mcrd
