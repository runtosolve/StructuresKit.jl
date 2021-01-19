module AISC360_16

export e3

function e3(Fe, Fy, Ag, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.80
    elseif ASDorLRFD==1
        StrengthFactor=0.85
    else
        StrengthFactor=1.0   #to just get nominal strength
    end

    if Fy/Fe <= 2.25
        Fcr = 0.658^(Fy/Fe) * Fy
    elseif λc > 2.25
        Fcr = 0.877 * Fe
    end

    Pn = Fcr * Ag

    ePn = Pn * StrengthFactor

    return Pn, ePn

end

end  #module