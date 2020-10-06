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


function h42(Mxbar,Mybar,Bbar,Mybar_freeflange, Maxℓo,Mayℓo,Ba, Mayℓo_freeflange)


    ActionMx = abs.(Mxbar./(Maxℓo))
    ActionMy = abs.(Mybar./(Mayℓo))
    ActionB = abs.(Bbar./(Ba))
    ActionMy_freeflange = abs.(Mybar_freeflange./(Mayℓo_freeflange))

    Interaction = ActionMx .+ ActionMy .+ ActionB .+ ActionMy_freeflange

    return ActionMx, ActionMy, ActionB, ActionMy_freeflange, Interaction

end


end #module
