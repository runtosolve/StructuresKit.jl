using Interpolations
using Polynomials


function ConcreteStressStrainModel(eps_co,eps_cu,f_c,f_t)

    #define concrete model strain range
    eps_c=range(0,eps_cu,length=100)

    #fill initial stress vector
    sig_c=zeros(Float64,length(eps_c))

    #define model for stress-strain relationship
    pre=findall(x->(x<=eps_co), eps_c)
    sig_c[pre]=f_c*(2*eps_c[pre]/eps_co-(eps_c[pre]./eps_co).^2)
    post=findall(x->(x>eps_co), eps_c)
    sig_c[post].=f_c.+((0.15.*f_c)./(eps_co.-eps_cu)).*(eps_c[post].-eps_co)

    #calculate initial elastic modulus
    E_ci=(sig_c[2]-sig_c[1])/(eps_c[2]-eps_c[1])

    #convert signs to represent compression
    eps_c=-eps_c
    sig_c=-sig_c

    #define concrete tension stress-strain curve
    f_cr=f_t  #tensile modulus of rupture
    eps_cr=f_cr/E_ci #strain associated with modulus of rupture

    eps_c=vcat(eps_cr,eps_c)
    # eps_c=vcat(eps_cr*1.000001,eps_c)
    # eps_c=vcat(0.003,eps_c)
    sig_c=vcat(f_cr,sig_c)
    # sig_c=vcat([0, 0],sig_c)

    return eps_c,sig_c,E_ci,eps_cr

end



function SteelStressStrainModel(E_si,sig_sy,sig_su,sig_sf,eps_sy,eps_sy1,eps_su,esp_sf)

    # define elastic part of steel stress-strain curve
    sig_s=[0 sig_sy]
    eps_s=[0 eps_sy]

    #define plateau and hardening part of steel stress-strain curve
    p = fit([eps_sy1, eps_su, eps_sf],[sig_sy, sig_su, sig_sf],2)

    sig_fit = zeros(11)
    eps_range = eps_sy1:(eps_sf-eps_sy1)./10:eps_sf

    for i=1:11  
    sig_fit[i] = p(eps_range[i])
    end

    eps_s=vcat(vec(eps_s),collect(range(eps_sy1,length=11,stop=eps_sf)))
    sig_s=vcat(vec(sig_s), sig_fit)

    # add compression part of steel stress-strain curve
    sig_s=vcat(-reverse(sig_s[2:end]),sig_s)
    eps_s=vcat(-reverse(eps_s[2:end]),eps_s)

    return eps_s, sig_s

end

function Response5000(M,P,h,s,w,A_s,z_s,E_ci,E_si,eps_cr,SteelStrain,SteelStress,ConcreteStrain,ConcreteStress)


    # define initial concrete locations in cross-section (from bottom, z=0)
    A_c=ones(s,1)*h/s*w
    z_c=collect(range(h/s*0.5,h-h/s*0.5,length=s))

    # initialize elastic modulus
    E_c=ones(s,1)*E_ci
    E_s=ones(length(A_s),1)*E_si

    SteelFit = LinearInterpolation(SteelStrain, SteelStress)
    ConcreteFit = LinearInterpolation(reverse(ConcreteStrain), reverse(ConcreteStress))

    #initialize stress, strain, and curvature
    eps_P=0
    eps_M_c=zeros(length(A_c),1)
    eps_M_s=zeros(length(A_s),1)
    rho_inv_hist=zeros(length(M))
    rho_inv=0

    sig_c_hist=zeros(length(A_c),1)
    sig_s_hist=zeros(length(z_s),1)
    eps_c_hist=zeros(length(A_c),1)
    eps_s_hist=zeros(length(z_s),1)

    for i=2:length(M)

        # E_c
        # E_s
        # eps_P
        # rho_inv

        # define load increment
        dP=P[i]-P[i-1]
        dM=M[i]-M[i-1]

        # calculate cross-section centroid
        z_bar=(sum(A_c.*z_c.*E_c)+sum(A_s.*z_s.*E_s))./((sum(A_c.*E_c)+sum(A_s.*E_s)))

        # calculate concrete and steel strains from P
        eps_P=eps_P+dP./(sum(E_s.*A_s)+sum(E_c.*A_c))

        # calculate 1 over the radius of curvature from M
        rho_inv=rho_inv.+dM./(sum(A_s.*E_s.*(z_s.-z_bar).^2).+sum(A_c.*E_c.*(z_c.-z_bar).^2))

        rho_inv_hist[i]=rho_inv

        # calculate strains from M
        eps_M_c=-(z_c.-z_bar).*rho_inv   #strain in concrete
        eps_M_s=-(z_s.-z_bar).*rho_inv   #strain in steel

        # combine axial and flexural strain, calculate strain history
        eps_c_hist=eps_M_c.+eps_P  #total strain in concrete
        eps_s_hist=eps_M_s.+eps_P  #total strain in steel

        index=findall(x->x<eps_cr,eps_c_hist)
        antindex=findall(x->x>=eps_cr,eps_c_hist)

        if abs(minimum(eps_c_hist))>abs(minimum(ConcreteStrain))
            println("The concrete has crushed, P=",P[i]," M=",M[i])
            return eps_c_hist,eps_s_hist,sig_c_hist,sig_s_hist,z_bar,z_c,rho_inv_hist
        end
        steelfailure=0
        for i=1:length(z_s)
            if abs(eps_s_hist[i])>maximum(SteelStrain)
                steelfailure=1
            else
                steelfailure=0
            end
        end

        if steelfailure==1
            println("The steel has reached its ultimate strain, P=",P[i], " M=",M[i])
            return eps_c_hist,eps_s_hist,sig_c_hist,sig_s_hist,z_bar,z_c,rho_inv_hist
        end

        #store previous stress and strain history
        sig_s_hist_old = sig_s_hist
        sig_c_hist_old = sig_c_hist
        eps_s_hist_old = eps_s_hist
        eps_c_hist_old = eps_c_hist

        # calculate concrete and steel stresses
        sig_c_hist[index]=ConcreteFit(eps_c_hist[index])
        if isempty(antindex)==false
             sig_c_hist[antindex].=0
        end
        sig_s_hist=SteelFit(eps_s_hist)


        #tangent modulus 
        E_c=(sig_c_hist .- sig_c_hist_old)./(eps_c_hist .- eps_c_hist_old)
        E_s=(sig_s_hist .- sig_s_hist_old)./(eps_s_hist .- eps_s_hist_old)
        # update elastic modulus
        # E_c=sig_c_hist./eps_c_hist
        # E_s=sig_s_hist./eps_s_hist

    end
    return eps_c_hist,eps_s_hist,sig_c_hist,sig_s_hist,z_bar,z_c,rho_inv_hist, E_c, E_s
end


# define concrete cross-section
w=50  #width
h=50  #depth/height
s=50*4  #number of cross-section divisions

# define rebar area and location (measured from bottom, z=0) in
# cross-section
     # number size location
     Rebar=[10 8 2.5]

RebarQuantity=Rebar[:,1]
RebarDiameter=Rebar[:,2]/8
RebarLocation=Rebar[:,3]

RebarSpacing=w./RebarQuantity

A_s=RebarQuantity.*pi.*(RebarDiameter./2).^2
z_s=RebarLocation

#define concrete compression engineering stress-strain model
#modified Hognestad
eps_co=0.002        #strain at f'c
eps_cu=0.0035       #ultimate strain at brittle failure
f_c=6.00            #concrete compressive strength
f_t=1/10*f_c        #concrete tensile strength

ConcreteStrain,ConcreteStress,E_ci,eps_cr=ConcreteStressStrainModel(eps_co,eps_cu,f_c,f_t)


E_si=29000                  #steel elastic modulus
sig_sy=60                   #steel yield stress
sig_su=80                   #steel ultimate stress
sig_sf=sig_sy*(1+0.10)      #steel fracture stress#
eps_sy=sig_sy/E_si          #steel yield strain
eps_sy1=sig_sy/E_si*(10)    #steel strain at end of yield plateau
eps_su=0.20                 #steel ultimate strain
eps_sf=0.30                 #steel fracture strain



SteelStrain,SteelStress=SteelStressStrainModel(E_si,sig_sy,sig_su,sig_sf,eps_sy,eps_sy1,eps_su,eps_sf)

M=range(0,38798.0,length=100)
P=range(0, 0,length=100)

eps_c_hist,eps_s_hist,sig_c_hist,sig_s_hist,z_bar,z_c,rho_inv_hist, E_c, E_s=Response5000(M,P,h,s,w,A_s,z_s,E_ci,E_si,eps_cr,SteelStrain,SteelStress,ConcreteStrain,ConcreteStress)

plot(eps_c_hist, z_c)

