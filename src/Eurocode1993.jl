module Eurocode1993

export eq13_64_2006, eq13_68_2006, eq13_618_2006, eq13_627_2006, eq13_628c_2006

function eq13_64_2006(Weff, fyb, γMO)

    McRd = Weff*fyb/γMO

end


function eq13_68_2006(hw, sw, ϕ, t, fyb, E, γMO)

    λw = 0.346 * sw/t * sqrt(fyb/E)

   if λw <= 0.83
        fbv = 0.58 * fyb
   elseif (λw > 0.83) & (λw <= 1.40)
        fbv = 0.48 * fyb / λw
   elseif λw > 1.40
        fbv = 0.67 * fyb / λw^2
   end

   VbRd = (hw/sin(deg2rad(ϕ)) *t * fbv) / γMO

   return VbRd, fbv, λw

end

function eq13_618_2006(α, t, fyb, E, r, la, ϕ, γM1)

    RwRd = α * t^2 * sqrt(fyb * E) * (1 - 0.1 * sqrt(r/t)) * (0.5 + sqrt(0.02 * la / t)) * (2.4 + (ϕ/90)^2) / γM1

end

function eq13_627_2006(NEd, NRd, MyEd, MyRd, MfRd, MplRd, VEd, VwRd)

     I_N = abs(NEd/NRd)
     I_M = abs(MyEd/MyRd)
     I_V = abs((1 - MfRd/MplRd)*(2*VEd/VwRd - 1)^2)
     I_total = I_N + I_M + I_V
    
     return I_N, I_M, I_V, I_total
end

function eq13_628c_2006(MEd, McRd, FEd, RwRd)

     I_M = abs(MEd/McRd)
     I_F = abs(FEd/RwRd)
     I_total = I_M + I_F

     return I_M, I_F, I_total

end

end #module




