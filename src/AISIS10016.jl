module AISIS10016


function e2(Fcre, Fy, Ag, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.80
    elseif ASDorLRFD==1
        StrengthFactor=0.85
    else
        StrengthFactor=1.0   #to just get nominal strength
    end

    λc = sqrt(Fy/Fcre)

    if λc <= 1.5
        Fn = 0.658^((λc)^2) * Fy
    else
        Fn = (0.877 /(λc)^2) * Fy
    end

    Pne = Ag * Fn

    ePne = Pne * StrengthFactor

    return Pn, ePne

end


function f321(Mne, Mcrℓ, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.67
    elseif ASDorLRFD==1
        StrengthFactor=0.90
    else
        StrengthFactor=1.0   #to just get nominal strength
    end

    λℓ=sqrt(Mne/Mcrℓ)

    if λℓ <= 0.776
        Mnℓ=Mne
    else
        Mnℓ=(1-0.15*(Mcrℓ/Mne)^0.4)*(Mcrℓ/Mne)^0.4*Mne
    end

    eMnℓ = Mnℓ * StrengthFactor

    return Mnℓ, eMnℓ

end

function f411(My, Mcrd, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.67
    elseif ASDorLRFD==1
        StrengthFactor=0.90
    else
        StrengthFactor=1.0   #to just get nominal strength
    end


    λd=sqrt(My/Mcrd)

    if λd <= 0.673
        Mnd=My
    else
        Mnd=(1-0.22*(Mcrd/My)^0.4)*(Mcrd/My)^0.4*My
    end

    eMnd=Mnd*StrengthFactor

    return Mnd, eMnd

end

function g21(E, h, t, Fy, Vcr, ASDorLRFD)

    if ASDorLRFD==0
        StrengthFactor=1/1.60
    elseif ASDorLRFD==1
        StrengthFactor=0.95
    else
        StrengthFactor=1.0   #to just get nominal strength
    end

    Aw = h.*t

    Vy=0.6 .*Aw .*Fy
    λv=sqrt.(Vy./Vcr)

    if λv <= 0.815
        Vn = Vy
    elseif (λv>0.815) & (λv<=1.227)
        Vn = 0.815 .*sqrt.(Vcr.*Vy)
    elseif λv > 1.227
        Vn = Vcr
    end

    eVn=Vn*StrengthFactor

    return Vn, eVn

end

function g231(h, t, Fcr)

    Aw = h.*t
    Vcr = Aw.*Fcr

    return Vcr

end

function g232(E, μ, kv, h, t)

    Fcr = (π^2 .*E.*kv)./(12 .*(1-μ.^2).*(h./t).^2)

end

function g233(a, h)

    if a./h <= 1.0
        kv = 4.00 .+ 5.34./(a./h).^2
    elseif a./h > 1.0
        kv = 5.34 .+ 4.00./(a./h).^2
    end

    return kv

end


function h21(Mbar, Vbar, Maℓo, Va)

    Interaction=sqrt.((Mbar./Maℓo).^2 .+ (Vbar./Va).^2)

end

function h121(Pbar, Mxbar, Mybar, Pa, Max, May)

    ActionP = abs.(Pbar./Pa)
    ActionMx = abs.(Mxbar./Max)
    ActionMy = abs.(Mybar./May)

    Interaction = ActionP .+ActionMx .+ ActionMy

    return ActionP, ActionMx, ActionMy, Interaction

end

function table23131(CorZ,t,b,d,θ)

    CorZ = convert(Int8, CorZ)

    Af=(b+d)*t
    Jf=1/3*b*t^3+1/3*d*t^3
    Cwf=0.0

    if CorZ==0

        Ixf=t*(t^2*b^2+4*b*d^3+t^2*b*d+t^2*b*d+d^4)/(12*(b+d))
        Iyf=t*(b^4+4*d*b^3)/(12*(b+d))
        Ixyf=(t*b^2*d^2)/(4*(b+d))
        xof=b^2/(2*(b+d))
        hxf=-(b^2+2*d*b)/(2*(b+d))
        hyf=-d^2/(2*(b+d))
        yof=hyf

    else

        Ixf=t*(t^2*b^2+4*b*d^3-4*b*d^3*cos(θ)^2+t^2*b*d+d^4-d^4*cos(θ)^2)/(12*(b+d))
        Iyf=t*(b^4+4*d*b^3+6*d^2*b^2*cos(θ)+4*d^3*b*cos(θ)^2+d^4*cos(θ)^2)/(12*(b+d))
        Ixyf=t*b*d^2*sin(θ)*(b+d*cos(θ))/(4*(b+d))
        xof=(b^2-d^2*cos(θ))/(2*(b+d))
        hxf=-(b^2+2*d*b+d^2*cos(θ))/(2*(b+d))
        hyf=-d^2*sin(θ)/(2*(b+d))
        yof=hyf
    end

    return Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof

end


function app23331(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, CurvatureSign, Lm, kϕ, Sf)

    θc=deg2rad(θc)
    bc=b-t/2-t/2*tan(θc/2)
    dc=d-t/2*tan(θc/2)

    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof=table23131(CorZ,t,bc,dc,θc)

    L = app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

    β=app23333(L, Lm, M1, M2, CurvatureSign)

    kϕfe=app23133(E,Ixf,xof,hxf,Cwf,Ixyf,Iyf,L,G,Jf)

    kϕwe=app23335(E,t,μ,ho,L)

    kϕfg=app23135(L,Af,xof,hxf,Ixyf,Iyf,yof,Ixf)

    kϕwg=app23336(f1,f2,ho,t,L)

    Fcrd=app23332(β, kϕfe, kϕwe, kϕ, kϕfg, kϕwg)

    Mcrd=Sf*Fcrd

    return Mcrd

end


function app23333(L, Lm, M1, M2)

    β=1+0.4(L/Lm)^0.7*(1+M1/M2)^0.7

    if β>=1.3
        β=1.3
    end

    if β<=1.0
        β=1.0
    end

    return β

end

function app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

    Lcrd=(4*π^4*ho*(1-μ^2)/t^3*(Ixf*(xof-hxf)^2+Cwf-Ixyf^2/Iyf*(xof-hxf)^2)+π^4*ho^4/720)^(1/4)

    L=minimum([Lcrd, Lm])

    return L

end

function app23133(E,Ixf,xof,hxf,Cwf,Ixyf,Iyf,L,G,Jf)

    kϕfe=(π/L)^4*(E*Ixf*(xof-hxf)^2+E*Cwf-E*Ixyf^2/Iyf*(xof-hxf)^2)+(π/L)^2*G*Jf

end

function app23335(E,t,μ,ho,L)

    kϕwe=E*t^3/(12*(1-μ^2))*(3/ho+(π/L)^2*19/60*ho+(π/L)^4*ho^3/240)

end

function app23135(L,Af,xof,hxf,Ixyf,Iyf,yof,Ixf)

    kϕfg=(π/L)^2*(Af*((xof-hxf)^2*(Ixyf/Iyf)^2-2*yof*(xof-hxf)*(Ixyf/Iyf)+hxf^2+yof^2)+Ixf+Iyf)

end

function app23336(f1,f2,ho,t,L)

    ξweb=(f1-f2)/f1

    kϕwg=ho*t*π^2/13440*((((45360*(1-ξweb)+62160)*(L/ho)^2)+448*π^2+(ho/L)^2*(53+3*(1-ξweb))*π^4)/(π^4+28*π^2*(L/ho)^2+420*(L/ho)^4))

    return kϕwg

end

function app23332(β, kϕfe, kϕwe, kϕ, kϕfg, kϕwg)

    Fcrd=β*(kϕfe+kϕwe+kϕ)/(kϕfg+kϕwg)

end

end #module
