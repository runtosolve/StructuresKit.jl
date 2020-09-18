module AISIS10024

using ..AISIS10016

export h411, h42


function h411(Cw, Fy, Wn, Bcrℓ, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.67
    elseif ASDorLRFD==1
        StrengthFactor=0.90
    else
        StrengthFactor=1.0   #to just get nominal strength
    end

    #calculate bimoment yield strength
    By=Cw.*Fy./Wn

    #calculate nominal bimoment strength including local buckling
    #for now use DSM local buckling flexural curve
    Bn, notused = AISIS10016.f321(By, Bcrℓ, ASDorLRFD)

    eBn = Bn * StrengthFactor

    return Bn, eBn

end


function h42(Mxbar,Mybar,Bbar,Maxℓo,Mayℓo,Ba)


    ActionMx = abs.(Mxbar./(Maxℓo))
    ActionMy = abs.(Mybar./(Mayℓo))
    ActionB = abs.(Bbar./(Ba))

    Interaction = ActionMx .+ ActionMy .+ ActionB

    return ActionMx, ActionMy, ActionB, Interaction

end


end #module
