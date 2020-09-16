using Plots
# using Plots.PlotMeasures
# using Calculus


#This program calculates the load-deformation response of a thin-walled column with discrete bracing.
#R0 May 17, 2019


#Solve as a boundary value problem in Julia.
using BoundaryValueDiffEq

#define all the variables for the thin-walled column
L=2438.0  #column length, mm
A=272  #column area, mm^2
xo=-32.59  #distance from centroid to shear center, x-direction, mm
yo=0       #distance from centroid to shear center, y-direction, mm
Ix=363370  #strong axis moment of inertia, mm^4
Iy=64100   #weak axis moment of inertia, mm^4
Io=Ix+Iy+(xo^2+yo^2)*A  #polar moment of inertia about the shear center, mm^4

E=200      #elastic modulus for steel, kN/mm^2
G=E/2.6    #shear modulus for steel, kN/mm^2

J=188   #St. Venant torsional constant, mm^4
Cw=122720891  #Warping torsion constant, mm^6

a1=L/1000  #weak axis flexural imperfection magnitude at midheight, mm
a2=L/1000  #strong axis flexural imperfection magnitude at midheight, mm
a3=0.00766 #twist imperfect magnitude at midheight, rad

hx=0   #location of x continuous translational spring from centroid, mm
hy=0   #location of y continous translational spring from centroid, mm

kx=0 #x continous translational spring stiffness, kN/m/m  (set high here to brace the weak axis)
ky=0       #y continous translational spring stiffness, kN/m/m
kϕ=200  #phi continous rotational spring stiffness, kN*m/rad/m

#define load range for column
Prange=0:5:15

# initialize midheight displacement vectors to be plotted against P
u_s_all=zeros(size(Prange))
v_s_all=zeros(size(Prange))
ϕ_s_all=zeros(size(Prange))

#define range of z along the column length
zspan = (0.0,L/2)

#define initial imperfection shape z ranges from 0 to 1
u_0(z) = a1*sin(pi*z/L)
v_0(z) = a2*sin(pi*z/L)
ϕ_0(z) = a3*sin(pi*z/L)

#write the system of ODEs for the column
function column!(dy,y,p,z)

    #define column deformations
    u  = y[1]   #deflection
    du = y[2]   #first derivative
    ddu = y[3]  #second derivative
    dddu= y[4]  #third derivative
    v  = y[5]   #deflection
    dv = y[6]   #first derivative
    ddv = y[7]  #second derivative
    dddv= y[8]  #third derivative
    ϕ = y[9]    #twist
    dϕ = y[10]  #first derivative
    ddϕ = y[11] #second derivative
    dddϕ = y[12] #third derivative

    # cdm1 = y[13]
    # cdm2 = y[14]
    # cdm3 = y[15]
    # cdm4 = y[16]
    # cdm5 = y[17]
    # cdm6 = y[18]

    u  = y[1]   #deflection

    #calculate second derivative of imperfections
    ddu_o=u_0''(z)
    ddv_o=v_0''(z)
    ddϕ_o=ϕ_0''(z)

    #define system of first order ODEs
    dy[1] = du
    dy[2] = ddu
    dy[3] = dddu
    dy[4] = (-P*(ddu+ddu_o))/(E*Iy)

    dy[5] = dv
    dy[6] = ddv
    dy[7] = dddv
    dy[8] = (-P*(ddv-xo*ddϕ+ddv_o-xo*ddϕ_o))/(E*Ix)

    dy[9] = dϕ
    dy[10] = ddϕ
    dy[11] = dddϕ
    dy[12] = ((G*J-P*Io/A)*ddϕ-P*yo*ddu+P*xo*ddv-P*yo*ddu+P*xo*ddv_o-P*(Io/A)*ddϕ_o)/(E*Cw)

#add dummy ODEs to allow more BC residuals
    # dy[13]=du
    # dy[14]=du
    # dy[15]=du
    # dy[16]=du
    # dy[17]=du
    # dy[18]=du
end

#define the boundary conditions for the column
function bc!(residual, y, p, z)
    residual[1] = y[1][1]       # u(0)=0
    # residual[2] = y[end][1]     # u(L)=0
    residual[2] = y[1][3]       # u''(0)=0   pinned end at x=0
    # residual[4] = y[end][3]     # u''(L)=0  pinned end at x=L
    residual[3] = y[1][5]       # v(0)=0
    # residual[6] = y[end][5]     # v(L)=0
    residual[4] = y[1][7]       # v''(0)=0   pinned end at x=0
    # residual[8] = y[end][7]     # v''(L)=0  pinned end at x=L
    residual[5] = y[1][9]       # ϕ(0)=0
    # residual[10] = y[end][9]    # ϕ(L)=0
    residual[6] = y[1][11]    # ϕ''(0)=0   pinned end at x=0
    # residual[12] = y[end][11]  # ϕ''(L)=0  pinned end at x=L

    # #transition conditions at midheight brace
    residual[7] = y[end][2]  #u'(L/2)=0
    residual[8] = y[end][6]  #v'(L/2)=0
    residual[9] = y[end][10] #ϕ'(L/2)=0
    residual[10] = E*Iy*y[end][4]-(kx*y[end][1]+kx*(yo-hy)*y[end][9])/2
    residual[11] = E*Ix*y[end][8]-(ky*y[end][5]-ky*(xo-hx)*y[end][9])/2
    residual[12] = E*Cw*y[end][12]-(kx*(yo-hy)*y[end][1]-ky*(xo-hx)*y[end][5]+kx*(yo-hy)^2*y[end][9]+ky*(xo-hx)^2*y[end][9]+kϕ*y[end][9])/2
end


bvp = BVProblem(column!, bc!, [a1,0,0,0,a2,0,0,0,a3,0,0,0], zspan)
sol = solve(bvp, GeneralMIRK4(), dt=0.01*L/2)


#loop over load range
for i in eachindex(Prange)

    global P=Prange[i]  #define compressive load on column

    #solve for the column displacements as a boundary value problem
    #define initial conditions for u, v, and ϕ as a1, a2, and a3
    bvp = BVProblem(column!, bc!, [a1,0,0,0,a2,0,0,0,a3,0,0,0], zspan)
    #dt= calculation increment distance along column, z/L
    sol = solve(bvp, GeneralMIRK4(), dt=0.01*L/2)

    #grab z range used by Julia
    global z=sol.t

    #grab column deformations as vectors
    global u_s=(s->s[1]).(sol.u)  #u
    global v_s=(s->s[5]).(sol.u)  #v
    global ϕ_s=(s->s[9]).(sol.u)  #ϕ

    #pull out displacements at column midheight
    u_s_all[i]=u_s[end]
    v_s_all[i]=v_s[end]
    ϕ_s_all[i]=ϕ_s[end]

end

#calculate total column deflection at midheight
u_total=u_s_all.+a1
v_total=v_s_all.+a2
ϕ_total=ϕ_s_all.+a3

#plot column displacements
plot(u_s.+a1,z,yaxis=("Column axial compression, P (kN)",font(12)),xaxis=("column deflection, u(z) (mm)",font(10)),legend=true)
plot!(v_s.+a2,z,yaxis=("distance along column, z (in.)"),xaxis=("column deflection, u(z) (mm)"),legend=true)
plot!(ϕ_s.+a3,z,yaxis=("distance along column, x (in.)",font(12)),xaxis=("column deflection, y(x) (in.)",font(12)),legend=true)

# #plot column displacements
# plot(u_total,Prange,yaxis=("Column axial compression, P (kN)",font(12)),xaxis=("column deflection, u(z) (mm)",font(10)),legend=true)
# plot!(v_total,Prange,yaxis=("distance along column, z (in.)"),xaxis=("column deflection, u(z) (mm)"),legend=true)
# plot!(ϕ_total,Prange,yaxis=("distance along column, x (in.)",font(12)),xaxis=("column deflection, y(x) (in.)",font(12)),legend=true)
