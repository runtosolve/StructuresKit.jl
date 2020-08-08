module AISIS10024

function h411(Cw, Fy, Wn, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.67
    elseif ASDorLRFD==1
        StrengthFactor=0.90
    else
        StrengthFactor=1.0   #to just get nominal strength
    end

    #calculate bimoment nominal strength
    Bn=Cw.*Fy./Wn*StrengthFactor

    return Bn

end


function h42(Mxbar,Mybar,Bbar,Maxℓo,Mayℓo,Ba)


    ActionMx = abs.(Mxbar./(Maxℓo))
    ActionMy = abs.(Mybar./(Mayℓo))
    ActionB = abs.(Bbar./(Ba))

    Interaction = ActionMx .+ ActionMy .+ ActionB

    return ActionMx, ActionMy, ActionB, Interaction

end


end #module
