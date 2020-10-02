module CrossSection

export CZcenterline, CUFSMtemplate, CUFSMproperties, CZflange_template, shear_stress, shear_force


function CZcenterline(H,B1,D1,q1,B2,D2,q2,ri1,ri2,ri3,ri4,t)

	#Translated from CUFSM v5.01 Matlab download, thanks Ben.

	#Calculate centerline dimensions of Cees and Zees from outside dimensions
	#and inner radii.

	# BWS 2015
	# reference AISI Design Manual for the lovely corner radius calcs.
	# For template calc, convert outer dimensons and inside radii to centerline
	# dimensiosn throughout
	# convert the inner radii to centerline if nonzero
	if ri1==0
	    r1=0
	else
	    r1=ri1+t/2
	end

	if ri2==0
	    r2=0
	else
	    r2=ri2+t/2
	end

	if ri3==0
	    r3=0
	else
	    r3=ri3+t/2
	end

	if ri4==0
	    r4=0
	else
	    r4=ri4+t/2
	end

	h=H-t/2-r1-r3-t/2

	if D1==0
	    b1=B1-r1-t/2
	    d1=0
	else
	    b1=B1-r1-t/2-(r2+t/2)*tan(q1/2)
	    d1=(D1-(r2+t/2)*tan(q1/2))
	end

	if D2==0
	    b2=B2-r3-t/2
	    d2=0
	else
	    b2=B2-r3-t/2-(r4+t/2)*tan(q2/2)
	    d2=(D2-(r4+t/2)*tan(q2/2))
	end

	return h, b1, d1, q1, b2, d2, q2, r1, r2, r3, r4, t

end


function CUFSMtemplate(CorZ,h,b1,b2,d1,d2,r1,r2,r3,r4,q1,q2,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)
	#Translated from CUFSM v5.01 Matlab download, thanks Ben.

	#BWS
	#August 23, 2000
	#2015 modification to allow for d1=d2=0 and creation of a track with same template
	#2015 addition to allow outer dimensions and inner radii to be used
	#2015 addition to control element discretization

	#CorZ=determines sign conventions for flange 1=C 2=Z
	if CorZ==2
		cz=-1;
	else
		cz=1;
	end
	#channel template

	#convert angles to radians
	q1=q1*π/180;
	q2=q2*π/180;
	#

	#if center is not 1 then outer dimensions and inner radii came in and these
	#need to be corrected to all centerline for the use of this template
	if center==1
	else
	    #label all the outer dimensions and inner radii
	    H=h;
	    B1=b1;
	    D1=d1;
	    ri1=r1;
	    ri2=r2;
	    B2=b2;
	    D2=d2;
	    ri3=r3;
	    ri4=r4;
	    h,b1,d1,q1,b2,d2,q2,r1,r2,r3,r4,t = CZcenterline(H,B1,D1,q1,B2,D2,q2,ri1,ri2,ri3,ri4,t);
	end


	#rest of the dimensions are "flat dimensions" and acceptable for modeling
	if (r1==0.0) & (r2==0.0) & (r3==0.0) & (r4==0.0)
	    if (d1==0.0) & (d2==0.0)
	        #track or unlipped Z with sharp corners
	        geom=[1 b1            0
	            2 0             0
	            3 0             h
	            4 cz*(b2)       h];
	        n=[nb1 nh nb2];
	    else
	        #lipped C or Z with sharp corners
	        geom=[1 b1+d1*cos(q1) d1*sin(q1)
	            2 b1            0
	            3 0             0
	            4 0             h
	            5 cz*(b2)       h
	            6 cz*(b2+d2*cos(q2)) h-d2*sin(q2)];
	        n=[nd1 nb1 nh nb2 nd2];
	    end
	else
	    if (d1==0.0) & (d2==0.0)
	        geom=[1 r1+b1                            0
	            2 r1                               0
	            3 0                                r1
	            4 0                                r1+h
	            5 cz*r3                               r1+h+r3
	            6 cz*(r3+b2)                            r1+h+r3];
	        n=[nb1 nr1 nh nr3 nb2];
	    else
	        geom=[1 r1+b1+r2*cos(π/2-q1)+d1*cos(q1) r2-r2*sin(π/2-q1)+d1*sin(q1)
	            2 r1+b1+r2*cos(π/2-q1)            r2-r2*sin(π/2-q1)
	            3 r1+b1                            0
	            4 r1                               0
	            5 0                                r1
	            6 0                                r1+h
	            7 cz*r3                               r1+h+r3
	            8 cz*(r3+b2)                            r1+h+r3
	            9 cz*(r3+b2+r4*cos(π/2-q2))            r1+h+r3-r4+r4*sin(π/2-q2)
	            10 cz*(r3+b2+r4*cos(π/2-q2)+d2*cos(q2)) r1+h+r3-r4+r4*sin(π/2-q2)-d2*sin(q2)];
	        n=[nd1 nr2 nb1 nr1 nh nr3 nb2 nr4 nd2];
	    end
	end
	#number of elements between the geom coordinates
	node=zeros(sum(n)+1,8)

	for i=1:size(geom)[1]-1

	    start=geom[i,2:3];
	    stop=geom[i+1,2:3];
	    if i==1
	        nstart=1;
	    else
	        nstart=sum(n[1:i-1])+1;
	    end

		node[nstart,:]=[nstart; geom[i,2:3]; 1; 1; 1; 1; 1.0];

	    if (r1==0.0) & (r2==0.0) & (r3==0.0) & (r4==0.0)
	        #------------------------
	        #SHARP CORNER MODEL
	        for j=1:n[i]-1
	            node[nstart+j,:]=[nstart+j; start+(stop-start)*j/n[i]; 1; 1; 1; 1; 1.0];
	        end
	        #------------------------
	    else
	        #ROUND CORNER MODEL
	        if (d1==0.0) & (d2==0.0) #track or unlipped Z geometry
	            #------------------------
	            #UNLIPPED C OR Z SECTION
	            if maximum(i.==[1 3 5]) #use linear interpolation
	                for j=1:n[i]-1
	                    node[nstart+j,:]=[nstart+j; start+(stop-start)*j/n[i]; 1; 1; 1; 1; 1.0];
	                end
	            else #we are in a corner and must be fancier
	                for j=1:n[i]-1
	                    if i==2
	                        r=r1;
							xc=r1;
							zc=r1;
							qstart=π/2;
							dq=π/2*j/n[i];
	                    end
	                    if i==4
	                        r=r3;
							xc=cz*r3;
							zc=r1+h;
							qstart=(1==cz)*π;
							dq=cz*π/2*j/n[i];
	                    end
	                    x2=xc+r*cos(qstart+dq);
	                    z2=zc-r*sin(qstart+dq); #note sign on 2nd term is negative due to z sign convention (down positive)
	                    node[nstart+j,:]=[nstart+j; x2; z2; 1; 1; 1; 1; 1.0];
	                end
	            end
	            #------------------------
	        else
	            #LIPPED C OR Z SECTION
	            #------------------------
	            if maximum(i.==[1 3 5 7 9]) #use linear interpolation
	                for j=1:n[i]-1
	                    node[nstart+j,:]=[nstart+j; start+(stop-start)*j/n[i]; 1; 1; 1; 1; 1.0];
	                end
	            else #we are in a corner and must be fancier
	                for j=1:n[i]-1
	                    if i==2
	                        r=r2;
							xc=r1+b1;
							zc=r2;
							qstart=π/2-q1;
							dq=q1*j/n[i];
	                    end
	                    if i==4
	                        r=r1;
							xc=r1;
							zc=r1;
							qstart=π/2;
							dq=π/2*j/n[i];
	                    end
	                    if i==6
	                        r=r3;
							xc=cz*r3;
							zc=r1+h;
							qstart=(1==cz)*π;
							dq=cz*π/2*j/n[i];
	                    end
	                    if i==8
	                        r=r4;
							xc=cz*(r3+b2);
							zc=r1+h+r3-r4;
							qstart=3*π/2;
							dq=cz*q2*j/n[i];
	                    end
	                    x2=xc+r*cos(qstart+dq);
	                    z2=zc-r*sin(qstart+dq); #note sign on 2nd term is negative due to z sign convention (down positive)
	                    node[nstart+j,:]=[nstart+j; x2; z2; 1; 1; 1; 1; 1.0];
	                end
	            end
	            #------------------------
	        end
	    end
	end



	#GET THE LAST NODE ASSIGNED
	if (r1==0.0) & (r2==0.0) & (r3==0.0) & (r4==0.0)
	    if (d1==0.0) & (d2==0.0)
	        i=4;
	        node[sum(n[1:i-1])+1,:]=[sum(n[1:i-1])+1; geom[i,2:3]; 1; 1; 1; 1; 1.0];
	    else
	        i=6;
	        node[sum(n[1:i-1])+1,:]=[sum(n[1:i-1])+1; geom[i,2:3]; 1; 1; 1; 1; 1.0];
	    end
	else
	    if (d1==0.0) & (d2==0.0)
	        i=6;
	        node[sum(n[1:i-1])+1,:]=[sum(n[1:i-1])+1; geom[i,2:3]; 1; 1; 1; 1; 1.0];
	    else
	        i=10;
	        node[sum(n[1:i-1])+1,:]=[sum(n[1:i-1])+1; geom[i,2:3]; 1; 1; 1; 1; 1.0];
	    end
	end

	elem=zeros(sum(n),5)
	for i=1:size(node)[1]-1
   		elem[i,:]=[i i i+1 t 100];
	end

	#set some default properties
	if kipin==1
		prop=[100 29500 29500 0.3 0.3 29500/(2*(1+0.3))];
	else
		prop=[100 203000 203000 0.3 0.3 203000/(2*(1+0.3))];
	end

	#set some default lengths
	big=maximum([h;b1;b2;d1;d2]);
	lengths=exp10.(range(log10(big/10), stop=log10(big*1000), length=50))

	springs=[0]
	constraints=[0]

	return prop,node,elem,lengths,springs,constraints,geom,cz

end





function CUFSMproperties(node, elem)

	#not matching CUFSM output exactly, need to check
	#because CUTWP function is used instead...

	#giving wrong answer for Cw and shear center when there is an isolated flange...

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

	return A, xc, yc, Ixc, Iyc, Ixyc, Imax, Imin, Th_p, Cw, J, Xs, Ys, w, Bx, By, B1, B2, VX, VY, Ai

end


#Get the node and element properties just for the flange and lip of a C or Z section.
#This code grabs the bottom flange and lip.
function CZflange_template(CorZ,H,Bc,Bt,Dc,Dt,r1,r2,r3,r4,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

    prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,Bc,Bt,Dc,Dt,r1,r2,r3,r4,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

    if CorZ == 2
        node[:, 2] = -node[:, 2]
    end

    index_xo = findall(x-> x==0.0, node[:, 2])
    index_yo = findall(y-> y==0.0, node[:, 3])
    index_o = intersect(index_xo, index_yo)
    index = 1:index_o[1]

	nodeflange = node[index,:]
    elemflange = elem[index[1:end-1],:]

    return nodeflange, elemflange

end

function shear_stress(V, Q, I, t)   #not working for Z section, only singly symmetric it seems

    τ = (V * Q) /(I * t)

    return τ

end

#integrate shear stress
function shear_force(τ, dA)

	V = 0.0
	for i = 1:length(τ) - 1
	    V = V + τ[i]*dA[i] + 1/2*(τ[i+1] - τ[i])*dA[i]
	end

	return V

end


end #module
