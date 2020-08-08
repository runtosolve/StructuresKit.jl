using StructuresKit

#test case
#8ZS2.25x070, no rotational spring, no moment gradient, Lm>>Lcrd
#sharp corners

CorZ = 1

t = 0.070

#out to out dimensions
ho = 8.0
b = 2.25
d = 0.930
θc = 50.0

E = 29500.0
μ = 0.30
G = E/(2*(1+μ))
f1 = 1.0
f2 = -1.0
M1 = 1.0
M2 = 1.0
Lm = 100.0
kϕ = 0.0

Sf = 9.4251/4.0

CurvatureSign = -1.0


Mcrd = AISIS10016.app23331(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, CurvatureSign, Lm, kϕ, Sf)

#AISI equations underpredict for a Z?
@test Mcrd ≈ 121.86 atol=0.09*Mcrd
