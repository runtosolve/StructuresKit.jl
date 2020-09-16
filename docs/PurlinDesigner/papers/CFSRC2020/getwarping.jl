using StructuresKit
using CSV, DataFrames


#read in cross-section dimensions from Gao and Moen (2013)
path = @__DIR__
filename = string(path, "/Gao_Moen_2013_Table1", ".csv")
dims = CSV.read(filename, DataFrame, header=false)



# for i in 1:49

i = 1

    #bring in cross-section dimensions
    Bc = dims[i,1]
    Dc = dims[i,2]
    θc = dims[i,4]

    Bt = dims[i,5]
    Dt = dims[i,6]
    θt = dims[i,8]

    H = dims[i,9]
    r = dims[i,10]
    t = dims[i,11]


    #calculate centerline dimensions

    CorZ = dims[i,13]+1
    center = 0  #outside dimensions
    kipin = 0

    nh = 4
    nb1 = 4
    nb2 = 4
    nd1 = 4
    nd2 = 4
    nr1 = 4
    nr2 = 4
    nr3 = 4
    nr4 = 4


    prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,Bc,Bt,Dc,Dt,r,r,r,r,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)


Cw, J, Xs, Ys, w, Bx, By, B1, B2 = warp(node, elem)


function warp(node, elem)

	#not matching CUFSM output exactly, need to check

    #Translated from CUFSM v5.01, thanks Ben.
	# %[A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22]=grosprop(node,elem)
	# %Program to compute the sectional properties of thin walled member.
	# %Cw, J , Warping function, Shear center location
	# % Badri Hiriyur Aug 20, 2002
	# %Input
	# %node=[# x z DOFX DOFZ DOFY DOF0 stress]
	# %elem=[# i j t]
	n=length(node[:,1])

	# %Basic Section Properties
	# % Area lengths, coordinates etc
	t=elem[:,4]
	x=node[:,2]
	y=node[:,3]

	A=0
	l=zeros(Float64, n-1)
	Ai=zeros(Float64, n-1)
	x_=zeros(Float64, n-1)
	y_=zeros(Float64, n-1)
	dx=zeros(Float64, n-1)
	dy=zeros(Float64, n-1)

	for i=1:n-1
	    l[i]=sqrt((x[i+1]-x[i])^2+(y[i+1]-y[i])^2)
	    Ai[i]=t[i]*l[i]
	    x_[i]=(1/2)*(x[i]+x[i+1])
	    y_[i]=(1/2)*(y[i]+y[i+1])
	    dx[i]=x[i+1]-x[i]
	    dy[i]=y[i+1]-y[i]
	    A=A+Ai[i]
	end

	# %Centroid
	xc=0
	yc=0
	for i=1:n-1
	    xc=xc+(1/A)*t[i]*l[i]*x_[i]
	    yc=yc+(1/A)*t[i]*l[i]*y_[i]
	end

	# %Moments of Inertia, Torsion constant
	J=0
	Ixc=0
	Iyc=0
	Ixyc=0
	for i=1:n-1
	    J=J+(1/3)*l[i]*t[i]^3
	    Ixc=Ixc+((y_[i]^2*Ai[i]+(1/12)*dy[i]^2*Ai[i]))
	    Iyc=Iyc+((x_[i]^2*Ai[i]+(1/12)*dx[i]^2*Ai[i]))
	    Ixyc=Ixyc+((x_[i]*y_[i]*Ai[i]+(1/12)*dx[i]*dy[i]*Ai[i]))
	end

	Ixc=Ixc-yc^2*A
	Iyc=Iyc-xc^2*A
	Ixyc=Ixyc-yc*xc*A
	# %Principal moments of inertia
	Imax=(1/2)*((Ixc+Iyc)+sqrt((Ixc-Iyc)^2+4*Ixyc^2))
	Imin=(1/2)*((Ixc+Iyc)-sqrt((Ixc-Iyc)^2+4*Ixyc^2))
	Th_p=1/2*(atan(-2*Ixyc,(Ixc-Iyc)))

	X = zeros(Float64, n)
	Y = zeros(Float64, n)
	# %Transform into new coordinates about principal axes
	for i=1:n
	    XY=[(x[i]-xc) (y[i]-yc); (y[i]-yc) -(x[i]-xc)]*[cos.(Th_p);sin.(Th_p)]
	    X[i]=XY[1]
	    Y[i]=XY[2]
	end

	# %Shearflow and Shear center
	VX=zeros(Float64, n)
	VY=zeros(Float64, n)
	dlX=zeros(Float64, n-1)
	dlY=zeros(Float64, n-1)
	dX = zeros(Float64, n-1)
	dY = zeros(Float64, n-1)

	for i=1:n-1
	    dX[i]=(1/l[i])*abs(X[i]*Y[i+1]-X[i+1]*Y[i])
	    dY[i]=(1/l[i])*abs(X[i]*Y[i+1]-X[i+1]*Y[i])

	    if (Y[i]*(X[i+1]-X[i]))<(X[i]*(Y[i+1]-Y[i]))
	        dlX[i]=1
	    elseif (Y[i]*(X[i+1]-X[i]))>(X[i]*(Y[i+1]-Y[i]))
	        dlX[i]=-1
	    elseif (Y[i]*(X[i+1]-X[i]))==(X[i]*(Y[i+1]-Y[i]))
	        dlX[i]=0
	    end

	    if (X[i]*(Y[i+1]-Y[i]))<(Y[i]*(X[i+1]-X[i]))
	        dlY[i]=1
	    elseif (X[i]*(Y[i+1]-Y[i]))>(Y[i]*(X[i+1]-X[i]))
	            dlY[i]=-1
	    elseif (X[i]*(Y[i+1]-Y[i]))==(Y[i]*(X[i+1]-X[i]))
	            dlY[i]=0
	    end

		VX[i+1]=VX[i]+Ai[i]*(Y[i+1]+Y[i])/2
	    VY[i+1]=VY[i]-Ai[i]*(X[i+1]+X[i])/2
	end

	Xs=0
	Ys=0
	if Imax!=0
	    for i=1:n-1
	        Xs=Xs+(-1/Imax)*dlX[i]*dX[i]*l[i]*(VX[i]+(1/6)*Ai[i]*(Y[i+1]+2*Y[i]))
	    end
	end
	if Imin!=0
	    for i=1:n-1
	        Ys=Ys+(-1/Imin)*dlY[i]*dY[i]*l[i]*(VY[i]-(1/6)*Ai[i]*(X[i+1]+2*X[i]))
	    end
	end

	# %Warping funcions and Warping Constant
	X_s = zeros(Float64, n)
	Y_s = zeros(Float64, n)
	ws = zeros(Float64, n)
	wa = zeros(Float64, n-1)
	ds = zeros(Float64, n-1)
	dls = zeros(Float64, n-1)

	X_s[1]=X[1]-Xs
	Y_s[1]=Y[1]-Ys
	ws[1]=0
	wa[1]=0
	ws_=0


	for i=1:n-1
	    X_s[i+1]=X[i+1]-Xs
	    Y_s[i+1]=Y[i+1]-Ys
	    ds[i]=(1/l[i])*abs(X_s[i]*Y_s[i+1]-X_s[i+1]*Y_s[i])
	    if (Y_s[i]*(X[i+1]-X[i]))<(X_s[i]*(Y[i+1]-Y[i]))
	        dls[i]=1
	    elseif (Y_s[i]*(X[i+1]-X[i]))>(X_s[i]*(Y[i+1]-Y[i]))
	        dls[i]=-1
	    elseif (Y_s[i]*(X[i+1]-X[i]))==(X_s[i]*(Y[i+1]-Y[i]))
	        dls[i]=0
	    end

	    ws[i+1]=ws[i]+ds[i]*l[i]*dls[i]
	    ws_=ws_+(1/A)*Ai[i]*(ws[i+1]+ws[i])/2
	end

	# xx[1]=0  not used it seems
	w = zeros(Float64, n)
	for i=1:n
	    w[i]=ws_-ws[i]
	end

	Cw=0
	for i=1:n-1
	    dw=w[i+1]-w[i]
	    wa=w[i]
	    Cw=Cw+t[i]*(wa^2*l[i]+(1/3)*dw^2*l[i]+wa*dw*l[i])
	end

	# %Monosymmetry Parameters Bx and By
	B1=0
	B2=0
	for i=1:n-1
	    Xa=X[i]
	    dX=X[i+1]-X[i]
	    Ya=Y[i]
	    dY=Y[i+1]-Y[i]
	    B1=B1+Ai[i]*1/12*(3*dY*dX^2+3*dY^3+4*Ya*dX^2+12*Ya*dY^2+8*dY*Xa*dX+12*Ya*Xa*dX+18*Ya^2*dY+6*dY*Xa^2+12*Ya*Xa^2+12*Ya^3)
	    B2=B2+Ai[i]*1/12*(3*dX^3+3*dX*dY^2+12*Xa*dX^2+4*Xa*dY^2+8*dX*Ya*dY+18*Xa^2*dX+12*Xa*Ya*dY+6*dX*Ya^2+12*Xa^3+12*Xa*Ya^2)
	end
	B1=(1/Imax)*B1-2*(Ys)
	B2=(1/Imin)*B2-2*(Xs)


	X=x .-xc.*ones(n)
	Y=y .-yc.*ones(n)
	Bx=0
	By=0
	for i=1:n-1
	    Xa=X[i]
	    dX=X[i+1]-X[i]
	    Ya=Y[i]
	    dY=Y[i+1]-Y[i]
	    Bx=Bx+Ai[i]*1/12*(3*dY*dX^2+3*dY^3+4*Ya*dX^2+12*Ya*dY^2+8*dY*Xa*dX+12*Ya*Xa*dX+18*Ya^2*dY+6*dY*Xa^2+12*Ya*Xa^2+12*Ya^3)
	    By=By+Ai[i]*1/12*(3*dX^3+3*dX*dY^2+12*Xa*dX^2+4*Xa*dY^2+8*dX*Ya*dY+18*Xa^2*dX+12*Xa*Ya*dY+6*dX*Ya^2+12*Xa^3+12*Xa*Ya^2)
	end
	Bx=(1/Ixc)*Bx-2*(Ys)
	By=(1/Iyc)*By-2*(Xs)

	return Cw, J, Xs, Ys, w, Bx, By, B1, B2

end





############






    testnum = string(i+1)


    filename = string("/Users/crismoen/RunToSolve/Papers/Moen_2020_CFSRC_PurlinDesigner", "/section_export", "/Test", testnum, "nodes", ".csv")
    CSV.write(filename,  DataFrame(node), header=false)

    filename = string("/Users/crismoen/RunToSolve/Papers/Moen_2020_CFSRC_PurlinDesigner", "/section_export", "/Test", testnum, "elements", ".csv")
    CSV.write(filename,  DataFrame(elem), header=false)

# end




using Plots
plot(node[:,2], node[:, 3], markershape = :circle)





kϕ = 1339  #N-mm/rad/mm   c=34 mm closest to test

#need
#section properties
#local buckling
