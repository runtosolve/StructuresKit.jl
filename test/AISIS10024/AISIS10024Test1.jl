using StructuresKit

#Test bimoment capacity equation against Bob Glauz's ballot example.
#see /testfiles/InternalForcesTest5/495A_Examples_v2.pdf
#check nominal and ASD calculations


Fy = 55
Cw = 11.9
Wn = 8.82

ASDorLRFD = 2   #nominal capacity
Bn = AISIS10024.h411(Cw, Fy, Wn, ASDorLRFD)

ASDorLRFD = 0   #ASD
BnASD = AISIS10024.h411(Cw, Fy, Wn, ASDorLRFD)


BnGlauz = 74.2
BnASDGlauz = 44.4

Error = abs((Bn - BnGlauz)/ BnGlauz)
ErrorASD = abs((BnASD - BnASDGlauz)/ BnASDGlauz)

#accept 1% error from numerical solution
@test Error <= 0.01
@test ErrorASD <= 0.01
