module CrossSection

using Statistics
using LinearAlgebra

export CUFSM_CZcenterline, CUFSMtemplate, CUFSMsection_properties, CZflange_template, shear_stress, shear_force



function CUFSM_CZcenterline(H,B1,D1,q1,B2,D2,q2,ri1,ri2,ri3,ri4,t)

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
	    h,b1,d1,q1,b2,d2,q2,r1,r2,r3,r4,t = CUFSM_CZcenterline(H,B1,D1,q1,B2,D2,q2,ri1,ri2,ri3,ri4,t);
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


function CUFSMsection_properties(coord,ends)

#Translated from Matlab to Julia on October 5, 2020.

#### Matlab notes below.

#Function modified for use in CUFSM by Ben Schafer in 2004 with permission
#of Sarawit. removed elastic buckling calcs & kept only section
#properties
#
#August 2005: additional modifications; program only handles
#singly-branched open sections; | single cell closed sections; arbitrary
#section designation added for other types.
#
#December 2006 bug fixes to B1 B2
#
#December 2015 extended to not crash on disconnected & arbitrary sections
#
#  Compute cross section properties
#----------------------------------------------------------------------------
#  Written by:
#       Andrew T. Sarawit
#       last revised:   Wed 10/25/01
#
#  Function purpose:
#       This function computes the cross section properties: area; centroid;
#       moment of inertia; torsional constant; shear center; warping constant;
#       B1; B2; elastic critical buckling load & the deformed buckling shape
#
#  Dictionary of Variables
#     Input Information:
#       coord[i,1:2]   .==  node i's coordinates
#                            coord[i,1] = X coordinate
#                            coord[i,2] = Y coordinate
#       ends[i,1:2]    .==  subelement i's nodal information
#                            ends[i,1] = start node #
#                            ends[i,2] = finish node #
#                            ends[i,3] = element's thickness
#       KL1            .==  effective unbraced length for bending about the 1-axis()
#       KL2            .==  effective unbraced length for bending about the 2-axis()
#       KL3            .==  effective unbraced length for twisting about the 3-axis()
#       force          .==  type of force applied
#                      .== "Pe"  : elastic critical axial force
#                      .== "Me1" : elastic critical moment about the 1-axis()
#                      .== "Me2" : elastic critical moment about the 2-axis()
#       exy[1:2]     .==  Pe eccentricities coordinates
#                            exy[1] = ex
#                            exy[2] = ey
#  Output Information:
#       A              .==  cross section area()
#       xc             .==  X coordinate of the centroid from orgin
#       yc             .==  Y coordinate of the centroid from origin
#       Ix             .==  moment of inertia about centroid X axes()
#       Iy             .==  moment of inertia about centroid Y axes()
#       Ixy            .==  product of inertia about centroid
#       Iz             .==  polar moment of inertia about centroid
#       theta          .==  rotation angle for the principal axes()
#       I1             .==  principal moment of inertia about centroid 1 axes()
#       I2             .==  principal moment of inertia about centroid 2 axes()
#       J              .==  torsional constant
#       xs             .==  X coordinate of the shear center from origin
#       ys             .==  Y coordinate of the shear center from origin
#       Cw             .==  warping constant
#       B1             .==  int[y*(x^2+y^2),s,0,L]   *BWS, x,y=prin. crd.
#       B2             .==  int[x*(x^2+y^2),s,0,L]
#                          where: x = x1+s/L*(x2-x1)
#                                 y = y1+s/L*(y2-y1)
#                                 L = lenght of the element
#       Pe[i]          .==  buckling mode i's elastic critical buckling load()
#       dcoord         .==  node i's coordinates of the deformed buckling shape
#                            coord[i,1,mode] = X coordinate
#                            coord[i,2,mode] = Y coordinate
#                          where: mode = buckling mode number
#
#  Note:
#     J;xs;ys;Cw;B1;B2;Pe;dcoord is not computed for close-section
#
#----------------------------------------------------------------------------
#
# find nele  .== total number of elements
#      nnode .== total number of nodes
#      j     .== total number of 2 element joints


    nele = size(ends,1);
    node = ends[:,1:2]; node = node[:]
    nnode = 0
    j = 0

    while isempty(node) == false
        i = findall(x-> x==node[1], node)
        deleteat!(node, i)
        # node[i] = []
        if size(i,1)==2
            j = j+1
        end
        nnode = nnode+1
    end

    # classify the section type()
    #This section modified in April 2006 by BWS to create an "arbitrary" category
    if j .== nele
        section = "close"; #single cell()
    elseif j .== nele-1
        section = "open";
    else
        section = "arbitrary"; #arbitrary section unidentified
            #in the future it would be good to handle multiple cross-sections in
            #one model; etc.; for now the code will bomb if more than a single()
            #section is used - due to inability to calculate section properties.
            #2015 decided to treat the secton as a fully composite section &
    end


    # if the section is close re-order the element
    if section == "close"
        xnele = (nele-1)
        for i = 1:xnele
            en = ends; en[i,2] = 0
            m,n = findall(ends[i,2]==en[:,1:2])
            if n==1
                ends[i+1,:] = en[m,:]
                ends[m,:] = en[i+1,:]
            elseif n .== 2
                ends[i+1,:] = en[m,[2 1 3]]
                ends[m,:] = en[i+1,[2 1 3]]
            end
        end
    end

    t = zeros(Float64, nele)
    xm = zeros(Float64, nele)
    ym = zeros(Float64, nele)
    xd = zeros(Float64, nele)
    yd = zeros(Float64, nele)
    L = zeros(Float64, nele)

    # find the element properties
    for i = 1:nele
        sn = ends[i,1]; fn = ends[i,2];

        sn = Int(sn)
        fn = Int(fn)

        # thickness of the element
        t[i] = ends[i,3]
        # compute the coordinate of the mid point of the element
        xm[i] = mean(coord[[sn fn],1])
        ym[i] = mean(coord[[sn fn],2])
        # compute the dimension of the element
        xd[i] = diff(vec(coord[[sn fn],1]))[1]
        yd[i] = diff(vec(coord[[sn fn],2]))[1]
        # compute the length of the element
        L[i] = norm([xd[i] yd[i]])
    end

    # compute the cross section area()
    A = sum(L.*t)
    # compute the centroid
    xc = sum(L.*t.*xm)/A
    yc = sum(L.*t.*ym)/A

    if abs(xc/sqrt(A)) .< 1e-12
        xc = 0
    end
    if abs(yc/sqrt(A)) .< 1e-12
        yc = 0
    end

    # compute the moment of inertia
    Ix = sum((yd.^2/12 .+(ym .-yc).^2).*L.*t)
    Iy = sum((xd.^2/12 .+(xm .-xc).^2).*L.*t)
    Ixy = sum((xd.*yd/12 .+(xm .-xc).*(ym .-yc)).*L.*t)

    if abs(Ixy/A^2) .< 1e-12
        Ixy = 0
    end

    # compute the rotation angle for the principal axes()
    theta = angle(Ix-Iy-2*Ixy*1im)/2

    coord12 = zeros(Float64, size(coord))

    # transfer the section coordinates to the centroid principal coordinates
    coord12[:,1] = coord[:,1] .-xc
    coord12[:,2] = coord[:,2] .-yc
    coord12 = [cos(theta) sin(theta); -sin(theta) cos(theta)]*transpose(coord12)
    coord12 = transpose(coord12)

    # find the element properties
    for i = 1:nele
        sn = ends[i,1]
        fn = ends[i,2]

        sn = Int(sn)
        fn = Int(fn)

        # compute the coordinate of the mid point of the element
        xm[i] = mean(coord12[[sn fn],1])
        ym[i] = mean(coord12[[sn fn],2])
        # compute the dimension of the element
        xd[i] = diff(vec(coord12[[sn fn],1]))[1]
        yd[i] = diff(vec(coord12[[sn fn],2]))[1]
    end

    # compute the principal moment of inertia
    I1 = sum((yd.^2/12+ym.^2).*L.*t)
    I2 = sum((xd.^2/12+xm.^2).*L.*t)

    if section == "close"
        # compute the torsional constant for close-section
        for i = 1:nele
            sn = ends[i,1]
            fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            p[i] = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
        end
        J = 4*sum(p.*L/2)^2/sum(L./t)
        xs = NaN; ys = NaN; Cw = NaN; B1 = NaN; B2 = NaN; Pe = NaN; dcoord = NaN; wn=NaN()
    elseif section == "open"
        # compute the torsional constant for open-section
        J = sum(L.*t.^3)/3

        # compute the shear center & initialize variables
        nnode = size(coord,1)
        w = zeros((nnode,2))
        w[Int(ends[1,1]),1] = ends[1,1]
        wo = zeros((nnode,2))
        wo[Int(ends[1,1]),1] = ends[1,1]
        Iwx = 0.0; Iwy = 0.0; wno = 0.0; Cw = 0.0

        for j = 1:nele
            i = 1
            while ((ends[i,1] in w[:,1]) & (ends[i,2] in w[:,1]))|(!(ends[i,1] in w[:,1])&(ends[i,2] ∉ w[:,1]))
                i = i+1
            end
            sn = ends[i,1]
            fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            p = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
            if w[sn,1]==0
                w[sn,1] = sn;
                w[sn,2] = w[fn,2]-p*L[i]
            elseif w[fn,1]==0
                w[fn,1] = fn;
                w[fn,2] = w[sn,2]+p*L[i]
            end
            Iwx = Iwx+(1/3*(w[sn,2]*(coord[sn,1]-xc)+w[fn,2]*(coord[fn,1]-xc))+1/6*(w[sn,2]*(coord[fn,1]-xc)+w[fn,2]*(coord[sn,1]-xc)))*t[i]* L[i];
            Iwy = Iwy+(1/3*(w[sn,2]*(coord[sn,2]-yc)+w[fn,2]*(coord[fn,2]-yc))+1/6*(w[sn,2]*(coord[fn,2]-yc)+w[fn,2]*(coord[sn,2]-yc)))*t[i]* L[i];
        end

        if (Ix*Iy-Ixy^2)!=0.0
            xs = (Iy*Iwy-Ixy*Iwx)/(Ix*Iy-Ixy^2)+xc
            ys = -(Ix*Iwx-Ixy*Iwy)/(Ix*Iy-Ixy^2)+yc
        else
            xs = xc; ys = yc
        end

        if abs(xs/sqrt(A)) .< 1e-12
            xs = 0
        end
        if abs(ys/sqrt(A)) .< 1e-12
            ys = 0
        end

        # compute the unit warping
        for j = 1:nele
            i = 1
            while ((ends[i,1] in wo[:,1]) & (ends[i,2] in wo[:,1]))|(!(ends[i,1] in wo[:,1])&(ends[i,2] ∉ wo[:,1]))
                i = i+1
            end
            sn = ends[i,1]
            fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            po = ((coord[sn,1]-xs)*(coord[fn,2]-ys)-(coord[fn,1]-xs)*(coord[sn,2]-ys))/L[i]
            if wo[sn,1]==0
                wo[sn,1] = sn;
                wo[sn,2] = wo[fn,2]-po*L[i]
            elseif wo[Int(ends[i,2]),1]==0
                wo[fn,1] = fn;
                wo[fn,2] = wo[sn,2]+po*L[i]
            end
            wno = wno+1/(2*A)*(wo[sn,2]+wo[fn,2])*t[i]* L[i];
        end
        wn = wno .-wo[:,2]

        # compute the warping constant
        for i = 1:nele
            sn = ends[i,1]; fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            Cw = Cw+1/3*(wn[sn]^2+wn[sn]*wn[fn]+wn[fn]^2)*t[i]* L[i]
        end

        # transfer the shear center coordinates to the centroid principal coordinates
        s12 = [cos(theta) sin(theta); -sin(theta) cos(theta)]* transpose([xs-xc ys-yc]);
        # compute the polar radius of gyration of cross section about shear center
        ro = sqrt((I1+I2)/A+s12[1]^2+s12[2]^2)

        # compute B1 & B2
        B1 = 0; B2 = B1
        for i = 1:nele
            sn = ends[i,1]
            fn = ends[i,2];

            sn = Int(sn)
            fn = Int(fn)

            x1 = coord12[sn,1]; y1 = coord12[sn,2]
            x2 = coord12[fn,1]; y2 = coord12[fn,2]
            B1 = B1+((y1+y2)*(y1^2+y2^2)/4+(y1*(2*x1^2+(x1+x2)^2)+y2*(2*x2^2+(x1+x2)^2))/12)*L[i]*t[i]
            B2 = B2+((x1+x2)*(x1^2+x2^2)/4+(x1*(2*y1^2+(y1+y2)^2)+x2*(2*y2^2+(y1+y2)^2))/12)*L[i]*t[i]
        end
        B1 = B1/I1-2*s12[2]
        B2 = B2/I2-2*s12[1];

        if abs(B1/sqrt(A)) .< 1e-12
            B1 = 0
        end
        if abs(B2/sqrt(A)) .< 1e-12
            B2 = 0
        end
    elseif section == "arbitrary"
        J = sum(L.*t.^3)/3
        xs = NaN; ys = NaN; Cw = NaN; B1 = NaN; B2 = NaN; Pe = NaN; dcoord = NaN; wn=NaN()

        #use the open section algorithm; modified to handle multiple parts; but
        #not completely generalized *that is a work for a future day* (17 Dec
        #2015) the primary goal here is to not crash the section property
        #calculators & applied stress generators when users have done a built
        #up section...

        # compute the torsional constant for open-section
        J = sum(L.*t.^3)/3

        # compute the shear center & initialize variables
        nnode = size(coord,1)
        w = zeros(nnode,2); w[Int(ends[1,1]),1] = ends[1,1]
        wo = zeros(nnode,2); wo[Int(ends[1,1]),1] = ends[1,1]
        Iwx = 0; Iwy = 0; wno = 0; Cw = 0

        for j = 1:nele
            i = 1
            while ((ends[i,1] in w[:,1]) & (ends[i,2] in w[:,1]))|(!(ends[i,1] in w[:,1])&(ends[i,2] ∉ w[:,1]))
                    i = i+1
                    if i>nele #brute force catch to continue calculation for multi part
                        i=nele
                        break
                    end
            end
            sn = ends[i,1]
            fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            p = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
            if w[sn,1]==0
                w[sn,1] = sn;
                w[sn,2] = w[fn,2]-p*L[i]
            elseif w[fn,1]==0
                w[fn,1] = fn;
                w[fn,2] = w[sn,2]+p*L[i]
            end
            Iwx = Iwx+(1/3*(w[sn,2]*(coord[sn,1]-xc)+w[fn,2]*(coord[fn,1]-xc))+1/6*(w[sn,2]*(coord[fn,1]-xc)+w[fn,2]*(coord[sn,1]-xc)))*t[i]* L[i];
            Iwy = Iwy+(1/3*(w[sn,2]*(coord[sn,2]-yc)+w[fn,2]*(coord[fn,2]-yc))+1/6*(w[sn,2]*(coord[fn,2]-yc)+w[fn,2]*(coord[sn,2]-yc)))*t[i]* L[i];
        end

        if (Ix*Iy-Ixy^2)~=0
            xs = (Iy*Iwy-Ixy*Iwx)/(Ix*Iy-Ixy^2)+xc
            ys = -(Ix*Iwx-Ixy*Iwy)/(Ix*Iy-Ixy^2)+yc
        else
            xs = xc; ys = yc
        end

        if abs(xs/sqrt(A)) .< 1e-12
            xs = 0
        end
        if abs(ys/sqrt(A)) .< 1e-12
            ys = 0
        end

        # compute the unit warping
        for j = 1:nele
            i = 1
            while ((ends[i,1] in wo[:,1]) & (ends[i,2] in wo[:,1]))|(!(ends[i,1] in wo[:,1])&(ends[i,2] ∉ wo[:,1]))
                    i = i+1
                    if i>nele #brute force catch to continue calculation for multi part
                        i=nele
                        break
                    end
            end
            sn = ends[i,1]; fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            po = ((coord[sn,1]-xs)*(coord[fn,2]-ys)-(coord[fn,1]-xs)*(coord[sn,2]-ys))/L[i]
            if wo[sn,1]==0
                wo[sn,1] = sn;
                wo[sn,2] = wo[fn,2]-po*L[i]
            elseif wo[ends[i,2],1]==0
                wo[fn,1] = fn;
                wo[fn,2] = wo[sn,2]+po*L[i]
            end
            wno = wno+1/(2*A)*(wo[sn,2]+wo[fn,2])*t[i]* L[i];
        end
        wn = wno-wo[:,2]

        # compute the warping constant
        for i = 1:nele
            sn = ends[i,1]; fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            Cw = Cw+1/3*(wn[sn]^2+wn[sn]*wn[fn]+wn[fn]^2)*t[i]* L[i]
        end

        # transfer the shear center coordinates to the centroid principal coordinates
        s12 = [cos(theta) sin(theta); -sin(theta) cos(theta)]*transpose([xs-xc ys-yc]);
        # compute the polar radius of gyration of cross section about shear center
        ro = sqrt((I1+I2)/A+s12[1]^2+s12[2]^2)

        # compute B1 & B2
        B1 = 0; B2 = B1
        for i = 1:nele
            sn = ends[i,1]; fn = ends[i,2];

            sn = Int(sn)
            fn = Int(fn)

            x1 = coord12[sn,1]; y1 = coord12[sn,2]
            x2 = coord12[fn,1]; y2 = coord12[fn,2]
            B1 = B1+((y1+y2)*(y1^2+y2^2)/4+(y1*(2*x1^2+(x1+x2)^2)+y2*(2*x2^2+(x1+x2)^2))/12)*L[i]*t[i]
            B2 = B2+((x1+x2)*(x1^2+x2^2)/4+(x1*(2*y1^2+(y1+y2)^2)+x2*(2*y2^2+(y1+y2)^2))/12)*L[i]*t[i]
        end
        B1 = B1/I1-2*s12[2]
        B2 = B2/I2-2*s12[1];

        if abs(B1/sqrt(A)) .< 1e-12
            B1 = 0
        end
        if abs(B2/sqrt(A)) .< 1e-12
            B2 = 0
        end

    end

    return A,xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,Cw,B1,B2,wn

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
