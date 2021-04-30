module CUFSM

using SparseArrays
using LinearAlgebra
using Statistics

export strip, stresgen, data, cutwp_prop2, templatecalc, template_out_to_in, SectionProperties

struct data

    prop::Matrix{Float64}
    node::Matrix{Float64}
    elem::Matrix{Float64}
    lengths::Array{Float64}
    springs::Matrix{Float64}
    constraints::Int64
    neigs::Int64
    curve::Vector{Matrix{Float64}}
    shapes::Vector{Matrix{Float64}}

end

struct SectionProperties

    node_geometry::Array{Float64, 2}
    element_info::Array{Float64, 2}
    A::Float64
    xc::Float64
    yc::Float64
    Ixx::Float64
    Iyy::Float64
    Ixy::Float64
    θ::Float64
    I1::Float64
    I2::Float64
    J::Float64
    xs::Float64
    ys::Float64
    Cw::Float64
    B1::Float64
    B2::Float64
    wn::Array{Float64, 1}

end



function template_out_to_in(H,B1,D1,q1,B2,D2,q2,ri1,ri2,ri3,ri4,t)

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


function templatecalc(CorZ,h,b1,b2,d1,d2,r1,r2,r3,r4,q1,q2,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)
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


function cutwp_prop2(coord,ends)

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
            en = deepcopy(ends); en[i,2] = 0
            index = findall(ends[i,2] .== en[:,1:2])  #update this for Julia, need to deal with CartesianIndex
            m=index[1][1]
            n=index[1][2]
            if n==1
                ends[i+1,:] = en[m,:]
                ends[m,:] = en[i+1,:]
            elseif n == 2
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

        p = zeros(Float64, nele)
        # compute the torsional constant for close-section
        for i = 1:nele
            sn = ends[i,1]
            fn = ends[i,2]

            sn = Int(sn)
            fn = Int(fn)

            p[i] = ((coord[sn,1]-xc)*(coord[fn,2]-yc)-(coord[fn,1]-xc)*(coord[sn,2]-yc))/L[i]
        end
        J = 4*sum(p.*L/2)^2/sum(L./t)
        xs = NaN; ys = NaN; Cw = NaN; B1 = NaN; B2 = NaN; Pe = NaN; dcoord = NaN; wn=[NaN; NaN]
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
        xs = NaN; ys = NaN; Cw = NaN; B1 = NaN; B2 = NaN; Pe = NaN; dcoord = NaN; wn=NaN

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

        if (Ix*Iy-Ixy^2) != 0
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
            elseif wo[Int(ends[i,2]),1]==0
                wo[fn,1] = fn;
                wo[fn,2] = wo[sn,2]+po*L[i]
            end
            wno = wno+1/(2*A)*(wo[sn,2]+wo[fn,2])*t[i]* L[i];
        end
        wn = wno .- wo[:,2]

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

    section_properties = SectionProperties(coord,ends,A,xc,yc,Ix,Iy,Ixy,theta,I1,I2,J,xs,ys,Cw,B1,B2,wn)

    return section_properties

end






function klocal(Ex,Ey,vx,vy,G,t,a,b,BC,m_a)
    #
    #Generate element stiffness matrix [k] in local coordinates

    # Inputs:
    # Ex;Ey;vx;vy;G: material properties
    # t: thickness of the strip [element]
    # a: length of the strip in longitudinal direction
    # b: width of the strip in transverse direction
    # BC: ["S-S"] a string specifying boundary conditions to be analyzed:
        #"S-S" simply-pimply supported boundary condition at loaded edges
        #"C-C" clamped-clamped boundary condition at loaded edges
        #"S-C" simply-clamped supported boundary condition at loaded edges
        #"C-F" clamped-free supported boundary condition at loaded edges
        #"C-G" clamped-guided supported boundary condition at loaded edges
    # m_a: longitudinal terms [or half-wave numbers] for this length()

    # Output:
    # k: local stiffness matrix; a totalm x totalm matrix of 8 by 8 submatrices.
    # k=[kmp]totalm x totalm block matrix
    # each kmp is the 8 x 8 submatrix in the DOF order [u1 v1 u2 v2 w1 theta1 w2 theta2]'

    # Z. Li June 2008
    # modified by Z. Li; Aug. 09; 2009
    # modified by Z. Li; June 2010

    E1=Ex/(1-vx*vy)
    E2=Ey/(1-vx*vy)
    Dx=Ex*t^3/(12*(1-vx*vy))
    Dy=Ey*t^3/(12*(1-vx*vy))
    D1=vx*Ey*t^3/(12*(1-vx*vy))
    Dxy=G*t^3/12
    #
    totalm = length(m_a); #Total number of longitudinal terms m
    #
    k=sparse(zeros(Float64, (8*totalm,8*totalm)))
    z0=zeros(4,4)
    for m=1:1:totalm
        for p=1:1:totalm
            #
            km_mp=zeros(4,4)
            kf_mp=zeros(4,4)
            um=m_a[m]*pi
            up=m_a[p]*pi
            c1=um/a
            c2=up/a
            #
            I1,I2,I3,I4,I5 = BC_I1_5(BC,m_a[m],m_a[p],a)
            #
            #asemble the matrix of Km_mp
            km_mp[1,1]=E1*I1/b+G*b*I5/3
            km_mp[1,2]=E2*vx*(-1/2/c2)*I3-G*I5/2/c2
            km_mp[1,3]=-E1*I1/b+G*b*I5/6
            km_mp[1,4]=E2*vx*(-1/2/c2)*I3+G*I5/2/c2

            km_mp[2,1]=E2*vx*(-1/2/c1)*I2-G*I5/2/c1
            km_mp[2,2]=E2*b*I4/3/c1/c2+G*I5/b/c1/c2
            km_mp[2,3]=E2*vx*(1/2/c1)*I2-G*I5/2/c1
            km_mp[2,4]=E2*b*I4/6/c1/c2-G*I5/b/c1/c2

            km_mp[3,1]=-E1*I1/b+G*b*I5/6
            km_mp[3,2]=E2*vx*(1/2/c2)*I3-G*I5/2/c2
            km_mp[3,3]=E1*I1/b+G*b*I5/3
            km_mp[3,4]=E2*vx*(1/2/c2)*I3+G*I5/2/c2

            km_mp[4,1]=E2*vx*(-1/2/c1)*I2+G*I5/2/c1
            km_mp[4,2]=E2*b*I4/6/c1/c2-G*I5/b/c1/c2
            km_mp[4,3]=E2*vx*(1/2/c1)*I2+G*I5/2/c1
            km_mp[4,4]=E2*b*I4/3/c1/c2+G*I5/b/c1/c2;
            km_mp=km_mp*t
            #
            #
            #asemble the matrix of Kf_mp
            kf_mp[1,1]=(5040*Dx*I1-504*b^2*D1*I2-504*b^2*D1*I3+156*b^4*Dy*I4+2016*b^2*Dxy*I5)/420/b^3
            kf_mp[1,2]=(2520*b*Dx*I1-462*b^3*D1*I2-42*b^3*D1*I3+22*b^5*Dy*I4+168*b^3*Dxy*I5)/420/b^3
            kf_mp[1,3]=(-5040*Dx*I1+504*b^2*D1*I2+504*b^2*D1*I3+54*b^4*Dy*I4-2016*b^2*Dxy*I5)/420/b^3
            kf_mp[1,4]=(2520*b*Dx*I1-42*b^3*D1*I2-42*b^3*D1*I3-13*b^5*Dy*I4+168*b^3*Dxy*I5)/420/b^3

            kf_mp[2,1]=(2520*b*Dx*I1-462*b^3*D1*I3-42*b^3*D1*I2+22*b^5*Dy*I4+168*b^3*Dxy*I5)/420/b^3
            kf_mp[2,2]=(1680*b^2*Dx*I1-56*b^4*D1*I2-56*b^4*D1*I3+4*b^6*Dy*I4+224*b^4*Dxy*I5)/420/b^3
            kf_mp[2,3]=(-2520*b*Dx*I1+42*b^3*D1*I2+42*b^3*D1*I3+13*b^5*Dy*I4-168*b^3*Dxy*I5)/420/b^3
            kf_mp[2,4]=(840*b^2*Dx*I1+14*b^4*D1*I2+14*b^4*D1*I3-3*b^6*Dy*I4-56*b^4*Dxy*I5)/420/b^3

            kf_mp[3,1]=kf_mp[1,3]
            kf_mp[3,2]=kf_mp[2,3]
            kf_mp[3,3]=(5040*Dx*I1-504*b^2*D1*I2-504*b^2*D1*I3+156*b^4*Dy*I4+2016*b^2*Dxy*I5)/420/b^3
            kf_mp[3,4]=(-2520*b*Dx*I1+462*b^3*D1*I2+42*b^3*D1*I3-22*b^5*Dy*I4-168*b^3*Dxy*I5)/420/b^3

            kf_mp[4,1]=kf_mp[1,4]
            kf_mp[4,2]=kf_mp[2,4]
            kf_mp[4,3]=(-2520*b*Dx*I1+462*b^3*D1*I3+42*b^3*D1*I2-22*b^5*Dy*I4-168*b^3*Dxy*I5)/420/b^3;#not symmetric
            kf_mp[4,4]=(1680*b^2*Dx*I1-56*b^4*D1*I2-56*b^4*D1*I3+4*b^6*Dy*I4+224*b^4*Dxy*I5)/420/b^3

            #assemble the membrane & flexural stiffness matrices
            kmp=[km_mp  z0
                   z0  kf_mp]
            #add it into local element stiffness matrix by corresponding to m
            k[8*(m-1)+1:8*m,8*(p-1)+1:8*p]=kmp
        end
    end

    return k

end


function BC_I1_5(BC,kk,nn,a)
    #
    # Calculate the 5 undetermined parameters I1;I2;I3;I4;I5 for local elastic
    # & geometric stiffness matrices.
    # BC: a string specifying boundary conditions to be analyzed:
        #"S-S" simply-pimply supported boundary condition at loaded edges
        #"C-C" clamped-clamped boundary condition at loaded edges
        #"S-C" simply-clamped supported boundary condition at loaded edges
        #"C-F" clamped-free supported boundary condition at loaded edges
        #"C-G" clamped-guided supported boundary condition at loaded edges
    #Outputs:
    #I1;I2;I3;I4;I5
        #calculation of I1 is the integration of Ym*Yn from 0 to a
        #calculation of I2 is the integration of Ym''*Yn from 0 to a
        #calculation of I3 is the integration of Ym*Yn'' from 0 to a
        #calculation of I3 is the integration of Ym*Yn'' from 0 to a
        #calculation of I4 is the integration of Ym"'*Yn'" from 0 to a
        #calculation of I5 is the integration of Ym"*Yn" from 0 to a

    if BC == "S-S"
        #For simply-pimply supported boundary condition at loaded edges
        if kk==nn
            I1=a/2
            I2=-kk^2*pi^2/a/2
            I3=-nn^2*pi^2/a/2
            I4=pi^4*kk^4/2/a^3
            I5=pi^2*kk^2/2/a
        else
            I1=0
            I2=0
            I3=0
            I4=0
            I5=0
        end
    elseif BC == "C-C"
        #For Clamped-clamped boundary condition at loaded edges
        #calculation of I1 is the integration of Ym*Yn from 0 to a
        if kk==nn
            if kk==1
                I1=3*a/8
            else()
                I1=a/4
            end
            I2=-(kk^2+1)*pi^2/4/a
            I3=-(nn^2+1)*pi^2/4/a
            I4=pi^4*((kk^2+1)^2+4*kk^2)/4/a^3
            I5=(1+kk^2)*pi^2/4/a
        else
            if kk-nn==2
                I1=-a/8
                I2=(kk^2+1)*pi^2/8/a-kk*pi^2/4/a
                I3=(nn^2+1)*pi^2/8/a+nn*pi^2/4/a
                I4=-(kk-1)^2*(nn+1)^2*pi^4/8/a^3
                I5=-(1+kk*nn)*pi^2/8/a
            elseif kk-nn==-2
                I1=-a/8
                I2=(kk^2+1)*pi^2/8/a+kk*pi^2/4/a
                I3=(nn^2+1)*pi^2/8/a-nn*pi^2/4/a
                I4=-(kk+1)^2*(nn-1)^2*pi^4/8/a^3
                I5=-(1+kk*nn)*pi^2/8/a
            else
                I1=0
                I2=0
                I3=0
                I4=0
                I5=0
            end
        end
    elseif (BC == "S-C") | (BC == "C-S")
        #For simply-clamped supported boundary condition at loaded edges
        #calculation of I1 is the integration of Ym*Yn from 0 to a
        if kk==nn
            I1=(1+(kk+1)^2/kk^2)*a/2
            I2=-(kk+1)^2*pi^2/a
            I3=-(kk+1)^2*pi^2/a
            I4=(kk+1)^2*pi^4*((kk+1)^2+kk^2)/2/a^3
            I5=(1+kk)^2*pi^2/a
        else
            if kk-nn==1
                I1=(kk+1)*a/2/kk
                I2=-(kk+1)*kk*pi^2/2/a
                I3=-(nn+1)^2*pi^2*(kk+1)/2/a/kk
                I4=(kk+1)*kk*(nn+1)^2*pi^4/2/a^3
                I5=(1+kk)*(1+nn)*pi^2/2/a
            elseif kk-nn==-1
                I1=(nn+1)*a/2/nn
                I2=-(kk+1)^2*pi^2*(nn+1)/2/a/nn
                I3=-(nn+1)*nn*pi^2/2/a
                I4=(kk+1)^2*nn*(nn+1)*pi^4/2/a^3
                I5=(1+kk)*(1+nn)*pi^2/2/a
            else()
                I1=0
                I2=0
                I3=0
                I4=0
                I5=0
            end
        end
        #
    elseif (BC == "C-F") | (BC = "F-C")
        #For clamped-free supported boundary condition at loaded edges
        #calculation of I1 is the integration of Ym*Yn from 0 to a
        if kk==nn
            I1=3*a/2-2*a*(-1)^(kk-1)/(kk-1/2)/pi
            I2=(kk-1/2)^2*pi^2*((-1)^(kk-1)/(kk-1/2)/pi-1/2)/a
            I3=(nn-1/2)^2*pi^2*((-1)^(nn-1)/(nn-1/2)/pi-1/2)/a
            I4=(kk-1/2)^4*pi^4/2/a^3
            I5=(kk-1/2)^2*pi^2/2/a
        else
            I1=a-a*(-1)^(kk-1)/(kk-1/2)/pi-a*(-1)^(nn-1)/(nn-1/2)/pi
            I2=(kk-1/2)^2*pi^2*((-1)^(kk-1)/(kk-1/2)/pi)/a
            I3=(nn-1/2)^2*pi^2*((-1)^(nn-1)/(nn-1/2)/pi)/a
            I4=0
            I5=0
        end
    elseif (BC == "C-G") | (BC == "G-C")
        #For clamped-guided supported boundary condition at loaded edges
        #calculation of I1 is the integration of Ym*Yn from 0 to a
        if kk==nn
            if kk==1
                I1=3*a/8
            else
                I1=a/4
            end
            I2=-((kk-1/2)^2+1/4)*pi^2/a/4
            I3=-((kk-1/2)^2+1/4)*pi^2/a/4
            I4=((kk-1/2)^2+1/4)^2*pi^4/4/a^3+(kk-1/2)^2*pi^4/4/a^3
            I5=(kk-1/2)^2*pi^2/a/4+pi^2/16/a
        else
            if kk-nn==1
                I1=-a/8
                I2=((kk-1/2)^2+1/4)*pi^2/a/8-(kk-1/2)*pi^2/a/8
                I3=((nn-1/2)^2+1/4)*pi^2/a/8+(nn-1/2)*pi^2/a/8
                I4=-nn^4*pi^4/8/a^3
                I5=-nn^2*pi^2/8/a
            elseif kk-nn==-1
                I1=-a/8
                I2=((kk-1/2)^2+1/4)*pi^2/a/8+(kk-1/2)*pi^2/a/8
                I3=((nn-1/2)^2+1/4)*pi^2/a/8-(nn-1/2)*pi^2/a/8
                I4=-kk^4*pi^4/8/a^3
                I5=-kk^2*pi^2/8/a
            else
                I1=0
                I2=0
                I3=0
                I4=0
                I5=0
            end
        end
    end

    return I1,I2,I3,I4,I5

end

function BC_I1_5_atpoint(BC,kk,nn,a,ys)
    #
    # Calculate the value of the longitundinal shape functions for discrete springs
    #
    # BC: a string specifying boundary conditions to be analyzed:
        #"S-S" simply-pimply supported boundary condition at loaded edges
        #"C-C" clamped-clamped boundary condition at loaded edges
        #"S-C" simply-clamped supported boundary condition at loaded edges
        #"C-F" clamped-free supported boundary condition at loaded edges
        #"C-G" clamped-guided supported boundary condition at loaded edges
    #Outputs:
    #I1;I5
        #calculation of I1 is the value of Ym[y/L]*Yn[y/L]
        #calculation of I5 is the value of Ym"(y/L)*Yn"(y/L)
    
    Ykk = Ym_at_ys(BC,kk,ys,a)
    Ynn = Ym_at_ys(BC,nn,ys,a)
    Ykkprime = Ymprime_at_ys(BC,kk,ys,a)    
    Ynnprime = Ymprime_at_ys(BC,nn,ys,a)    
    I1=Ykk*Ynn
    I5=Ykkprime*Ynnprime
    
    return I1, I5
            
end

function Ym_at_ys(BC,m,ys,a)
    #Longitudinal shape function values
    #written by BWS in 2015
    #could be called in lots of places, but now (2015) is hardcoded by Zhanjie
    #in several places in the interface
    #written in 2015 because wanted it for a new idea on discrete springs
    
    if BC == "S-S"
        Ym=sin(m*pi*ys/a)
    elseif BC == "C-C"
        Ym=sin(m*pi*ys/a) .*sin(pi*ys/a)
    elseif (BC == "S-C") | (BC == "C-S")
        Ym=sin((m+1)*pi*ys/a)+(m+1)/m*sin(m*pi*ys/a)
    elseif (BC == "C-F") | (BC == "F-C")
        Ym=1-cos((m-0.5)*pi*ys/a)
    elseif (BC == "C-G") | (BC = "G-C")
        Ym=sin((m-0.5)*pi*ys/a) .*sin(pi*ys/a/2)
    end
   
    return Ym

end


function Ymprime_at_ys(BC,m,ys,a)
    #First Derivative of Longitudinal shape function values
    #written by BWS in 2015
    #could be called in lots of places, but now (2015) is hardcoded by Zhanjie
    #in several places in the interface
    #written in 2015 because wanted it for a new idea on discrete springs
    
    if BC == "S-S"
        Ymprime=(pi*m*cos((pi*m*ys)/a))/a
    elseif BC == "C-C"
        Ymprime=(pi*cos((pi*ys)/a)*sin((pi*m*ys)/a))/a + (pi*m*sin((pi*ys)/a)*cos((pi*m*ys)/a))/a
    elseif (BC == "S-C") | (BC == "C-S")
        Ymprime=(pi*cos((pi*ys*(m + 1))/a)*(m + 1))/a + (pi*cos((pi*m*ys)/a)*(m + 1))/a
    elseif (BC == "C-F") | (BC == "F-C")
        Ymprime=(pi*sin((pi*ys*(m - 1/2))/a)*(m - 1/2))/a
    elseif (BC = "C-G") | (BC == "G-C")
        Ymprime=(pi*sin((pi*ys*(m - 1/2))/a)*cos((pi*ys)/(2*a)))/(2*a) + (pi*cos((pi*ys*(m - 1/2))/a)*sin((pi*ys)/(2*a))*(m - 1/2))/a
    end
    
    return Ymprime
            
end


function kglocal(a,b,Ty1,Ty2,BC,m_a)
    #
    #Generate geometric stiffness matrix [kg] in local coordinates

    # Inputs:
    # a: length of the strip in longitudinal direction
    # b: width of the strip in transverse direction
    # Ty1; Ty2: node stresses
    # BC: a string specifying boundary conditions to be analyzed:
        #"S-S" simply-pimply supported boundary condition at loaded edges
        #"C-C" clamped-clamped boundary condition at loaded edges
        #"S-C" simply-clamped supported boundary condition at loaded edges
        #"C-F" clamped-free supported boundary condition at loaded edges
        #"C-G" clamped-guided supported boundary condition at loaded edges
    # m_a: longitudinal terms [or half-wave numbers] for this length()

    # Output:
    # kg: local geometric stiffness matrix; a totalm x totalm matrix of 8 by 8 submatrices.
    # kg=[kgmp]totalm x totalm block matrix
    # each kgmp is the 8 x 8 submatrix in the DOF order [u1 v1 u2 v2 w1 theta1
    # w2 theta2]'

    # Z. Li; June 2008
    # modified by Z. Li; Aug. 09; 2009
    # modified by Z. Li; June 2010

    totalm = length(m_a); #Total number of longitudinal terms m
    kg=sparse(zeros(Float64, (8*totalm,8*totalm)))
    #
    for m=1:1:totalm
        for p=1:1:totalm
            #
            gm_mp=zeros(4,4)
            z0=zeros(4,4)
            gf_mp=zeros(4,4)
            um=m_a[m]*pi
            up=m_a[p]*pi
            #
            I1,I2,I3,I4,I5 = BC_I1_5(BC,m_a[m],m_a[p],a)
            #
            #asemble the matrix of gm_mp [symmetric membrane stability matrix]
            gm_mp[1,1]=b*(3*Ty1+Ty2)*I5/12
            gm_mp[1,3]=b*(Ty1+Ty2)*I5/12
            gm_mp[3,1]=gm_mp[1,3]
            gm_mp[2,2]=b*a^2*(3*Ty1+Ty2)*I4/12/um/up
            gm_mp[2,4]=b*a^2*(Ty1+Ty2)*I4/12/um/up
            gm_mp[4,2]=gm_mp[2,4]
            gm_mp[3,3]=b*(Ty1+3*Ty2)*I5/12
            gm_mp[4,4]=b*a^2*(Ty1+3*Ty2)*I4/12/um/up
            #
            #asemble the matrix of gf_mp [symmetric flexural stability matrix]
            gf_mp[1,1]=(10*Ty1+3*Ty2)*b*I5/35
            gf_mp[1,2]=(15*Ty1+7*Ty2)*b^2*I5/210/2
            gf_mp[2,1]=gf_mp[1,2]
            gf_mp[1,3]=9*(Ty1+Ty2)*b*I5/140
            gf_mp[3,1]=gf_mp[1,3]
            gf_mp[1,4]=-(7*Ty1+6*Ty2)*b^2*I5/420
            gf_mp[4,1]=gf_mp[1,4]
            gf_mp[2,2]=(5*Ty1+3*Ty2)*b^3*I5/2/420
            gf_mp[2,3]=(6*Ty1+7*Ty2)*b^2*I5/420
            gf_mp[3,2]=gf_mp[2,3]
            gf_mp[2,4]=-(Ty1+Ty2)*b^3*I5/140/2
            gf_mp[4,2]=gf_mp[2,4]
            gf_mp[3,3]=(3*Ty1+10*Ty2)*b*I5/35
            gf_mp[3,4]=-(7*Ty1+15*Ty2)*b^2*I5/420
            gf_mp[4,3]=gf_mp[3,4]
            gf_mp[4,4]=(3*Ty1+5*Ty2)*b^3*I5/420/2
            #assemble the membrane & flexural stiffness matrices
            kgmp=[gm_mp  z0
                   z0  gf_mp]
            #add it into local geometric stiffness matrix by corresponding to m
            kg[8*(m-1)+1:8*m,8*(p-1)+1:8*p]=kgmp;
        end
    end

    return kg

end

function elemprop(node,elem,nnodes,nelems)
    #BWS
    #1998 [last modified]
    #
    #INPUT
    #node: [node# x z dofx dofz dofy dofrot stress] nnodes x 8
    #elem: [elem# nodei nodej t] nelems x 4
    #OUTPUT
    #elprop: [elem# width alpha]
    #
    elprop=zeros(nelems,3)
    #
    for i=1:nelems
    	nodei = Int(elem[i,2])
    	nodej = Int(elem[i,3])
    	   xi = node[nodei,2]
    	   zi = node[nodei,3]
    	   xj = node[nodej,2]
    	   zj = node[nodej,3]
    	   dx = xj - xi
    	   dz = zj - zi
    	width = sqrt(dx^2 + dz^2)
    	alpha = atan(dz,dx)
    	elprop[i,:]=[i width alpha]
    end

    return elprop
end


function assemble(K,Kg,k,kg,nodei,nodej,nnodes,m_a)
    #
    #Add the element contribution to the global stiffness matrix

    #Outputs:
    # K: global elastic stiffness matrix
    # Kg: global geometric stiffness matrix
    # K & Kg: totalm x totalm submatrices. Each submatrix is similar to the
    # one used in original CUFSM for single longitudinal term m in the DOF order
    #[u1 v1...un vn w1 01...wn 0n]m'.

    # Z. Li; June 2008
    # modified by Z. Li; Aug. 09; 2009
    # Z. Li; June 2010

    totalm = length(m_a); #Total number of longitudinal terms m
    K2=sparse(zeros(Float64, (4*nnodes*totalm,4*nnodes*totalm)))
    K3=sparse(zeros(Float64, (4*nnodes*totalm,4*nnodes*totalm)))
    skip=2*nnodes
    for i=1:1:totalm
        for j=1:1:totalm
            #Submatrices for the initial stiffness
            k11=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+1:8*(j-1)+2]
            k12=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+3:8*(j-1)+4]
            k13=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+5:8*(j-1)+6]
            k14=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+7:8*(j-1)+8]
            k21=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+1:8*(j-1)+2]
            k22=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+3:8*(j-1)+4]
            k23=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+5:8*(j-1)+6]
            k24=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+7:8*(j-1)+8]
            k31=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+1:8*(j-1)+2]
            k32=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+3:8*(j-1)+4]
            k33=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+5:8*(j-1)+6]
            k34=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+7:8*(j-1)+8]
            k41=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+1:8*(j-1)+2]
            k42=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+3:8*(j-1)+4]
            k43=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+5:8*(j-1)+6]
            k44=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+7:8*(j-1)+8]
            #
            K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k11
            K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k12
            K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k21
            K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k22
            #
            K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k33
            K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k34
            K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k43
            K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k44
            #
            K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k13
            K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k14
            K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k23
            K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k24
            #
            K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k31
            K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k32
            K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k41
            K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k42
            #
            #Submatrices for the initial stiffness
            kg11=kg[8*(i-1)+1:8*(i-1)+2,8*(j-1)+1:8*(j-1)+2]
            kg12=kg[8*(i-1)+1:8*(i-1)+2,8*(j-1)+3:8*(j-1)+4]
            kg13=kg[8*(i-1)+1:8*(i-1)+2,8*(j-1)+5:8*(j-1)+6]
            kg14=kg[8*(i-1)+1:8*(i-1)+2,8*(j-1)+7:8*(j-1)+8]
            kg21=kg[8*(i-1)+3:8*(i-1)+4,8*(j-1)+1:8*(j-1)+2]
            kg22=kg[8*(i-1)+3:8*(i-1)+4,8*(j-1)+3:8*(j-1)+4]
            kg23=kg[8*(i-1)+3:8*(i-1)+4,8*(j-1)+5:8*(j-1)+6]
            kg24=kg[8*(i-1)+3:8*(i-1)+4,8*(j-1)+7:8*(j-1)+8]
            kg31=kg[8*(i-1)+5:8*(i-1)+6,8*(j-1)+1:8*(j-1)+2]
            kg32=kg[8*(i-1)+5:8*(i-1)+6,8*(j-1)+3:8*(j-1)+4]
            kg33=kg[8*(i-1)+5:8*(i-1)+6,8*(j-1)+5:8*(j-1)+6]
            kg34=kg[8*(i-1)+5:8*(i-1)+6,8*(j-1)+7:8*(j-1)+8]
            kg41=kg[8*(i-1)+7:8*(i-1)+8,8*(j-1)+1:8*(j-1)+2]
            kg42=kg[8*(i-1)+7:8*(i-1)+8,8*(j-1)+3:8*(j-1)+4]
            kg43=kg[8*(i-1)+7:8*(i-1)+8,8*(j-1)+5:8*(j-1)+6]
            kg44=kg[8*(i-1)+7:8*(i-1)+8,8*(j-1)+7:8*(j-1)+8]
            #
            K3[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=kg11
            K3[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=kg12
            K3[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=kg21
            K3[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=kg22
            #
            K3[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=kg33
            K3[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=kg34
            K3[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=kg43
            K3[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=kg44
            #
            K3[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=kg13
            K3[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=kg14
            K3[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=kg23
            K3[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=kg24
            #
            K3[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=kg31
            K3[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=kg32
            K3[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=kg41
            K3[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=kg42
        end
    end
    K=K+K2
    Kg=Kg+K3

    return K, Kg

end

function addspring(K,springs,nnodes,a,BC,m_a)
    #BWS
    #August 2000

    #[K]=addspring(K,springs,nnodes,a,m_a,BC)
    #Add spring stiffness to global elastic stiffness matrix
    #
    #K is the complete elastic stiffness matrix
    #springs is the definiton of any external springs added to the member
    #springs=[node# DOF[x=1,z=2,y=3,theta=4] ks]

    # modified by Z. Li; Aug. 09; 2009 for general B.C.
    # Z. Li; June 2010

    if (springs==0) | (isempty(springs))
        #nothing to calculate
    else
        totalm = length(m_a); #Total number of longitudinal terms m

        for i=1:size(springs,1)
            node=springs[i,1]
            dof=springs[i,2]
            k=springs[i,3]
            kflag=springs[i,4]

            if dof==1
                rc=2*node-1
            elseif dof==2
                rc=2*nnodes+2*node-1
            elseif dof==3
                rc=2*node
            elseif dof==4
                rc=2*nnodes+2*node
            else
                rc=1
                ks=0
            end
            for kk=1:1:totalm
                for nn=1:1:totalm
                    if kflag==0
                        ks=k; #k is the total stiffness & may be added directly
                    else
                        if dof==3
                            ks=0; #axial dof with a foundation stiffness has no net stiffness
                        else
                            I1,I2,I3,I4,I5 = BC_I1_5(BC,m_a[kk],m_a[nn],a)
                            ks=k*I1; #k is a foundation stiffness & an equivalent total stiffness must be calculated
                        end
                    end

                    K[4*nnodes*(kk-1)+rc,4*nnodes*(nn-1)+rc]=K[4*nnodes*(kk-1)+rc,4*nnodes*(nn-1)+rc]+ks
                end
            end
        end
    end

    return K

end

function trans(alpha,k,kg,m_a)
    #
    # Transfer the local stiffness into global stiffness
    # Zhanjie 2008
    # modified by Z. Li; Aug. 09; 2009

    totalm = length(m_a); #Total number of longitudinal terms m
    a=alpha

    #added by CDM
    num_dof = size(k)[1]
    gamma = zeros(Float64, (num_dof, num_dof))

    #
    z0=0
    gam=[cos(a)   z0    z0	  z0	-sin(a)	z0   z0	   z0
    	    z0	  1	    z0	  z0	  z0    z0	 z0    z0
    	    z0	  z0  cos(a)  z0	  z0	z0 -sin(a) z0
    	    z0	  z0	z0	   1	  z0	z0	 z0	   z0
    	   sin(a) z0	z0	  z0	cos(a)	z0	 z0	   z0
    	    z0	  z0	z0	  z0	  z0	1	 z0    z0
    	    z0	  z0  sin(a)  z0	  z0	z0	cos(a) z0
    	    z0	  z0	z0	  z0	  z0	z0	 z0	   1 ]
    #extend to multi-m
    for i=1:totalm
        gamma[8*(i-1)+1:8*i,8*(i-1)+1:8*i] = gam
    end
    #
    kglobal=gamma*k*gamma'
    kgglobal=gamma*kg*gamma'

    return kglobal,kgglobal
end


function constr_BCFlag(node,constraints)
    #
    # this subroutine is to determine flags for user constraints & internal [at node] B.transpose(C)s

    #inputs:
    #node: [node# x z dofx dofz dofy dofrot stress] nnodes x 8
    #constraints:: [node#e dofe coeff node#k dofk] e=dof to be eliminated k=kept dof dofe_node = coeff*dofk_nodek

    #Output:
    #BCFlag: 1 if there are user constraints | node fixities
           # 0 if there is no user constraints & node fixities

    # Z. Li; June 2010

    #Check for boundary conditions on the nodes
    nnodes = length(node[:,1])
    BCFlag=0
    for i=1:nnodes
        for j=4:7
            if node[i,j]==0
                BCFlag=1
                return
            end
        end
    end
    #Check for user defined constraints too
    if (isempty(constraints) | constraints==0) & (BCFlag==0)
        BCFlag=0
    else
        BCFlag=1
    end

    return BCFlag

end



function constr_user(node,cnstr,m_a)
    #
    #this routine creates the constraint matrix; Ruser; as defined by the user
    #
    #
    #input/output data
    #   node - same as elsewhere throughout this program
    #   cnstr - same as "constraints" throughout this program
    #   m_a - longitudinal terms to be included for this length()

    #   Ruser - the constraint matrix [in other words: base vectors] so that
    #               displ_orig = Ruser * displ_new

    # S. Adany; Feb 26; 2004
    # Z. Li; Aug 18; 2009 for general b.c.
    # Z. Li; June 2010

    nnode=length(node[:,1])
    ndof_m=4*nnode
    DOFreg=zeros(ndof_m,1)+1
    totalm = length(m_a); #Total number of longitudinal terms m
    # Ruser=I
    for ml=1:totalm
        #
        Ruser_m=I
        #to consider free DOFs
        for i=1:nnode
            for j=4:7
                if node[i,j]==0
                    if j==4
                        dofe=(i-1)*2+1
                    end
                    if j==6
                        dofe=i*2
                    end
                    if j==5
                        dofe=nnode*2+(i-1)*2+1
                    end
                    if j==7
                        dofe=nnode*2+i*2
                    end
                    DOFreg[dofe,1]=0
                end
            end
        end
        #
        #to consider master-slave constraints
        for i=1:length(cnstr[:,1])
            if length(cnstr[i,:])>=5
                #
                #nr of eliminated DOF
                nodee=cnstr[i,1]
                if cnstr[i,2]==1
                    dofe=(nodee-1)*2+1
                end
                if cnstr[i,2]==3
                    dofe=nodee*2
                end
                if cnstr[i,2]==2
                    dofe=nnode*2+(nodee-1)*2+1
                end
                if cnstr[i,2]==4
                    dofe=nnode*2+nodee*2
                end
                #
                #nr of kept DOF
                nodek=cnstr[i,4]
                if cnstr[i,5]==1
                    dofk=(nodek-1)*2+1
                end
                if cnstr[i,5]==3
                    dofk=nodek*2
                end
                if cnstr[i,5]==2
                    dofk=nnode*2+(nodek-1)*2+1
                end
                if cnstr[i,5]==4
                    dofk=nnode*2+nodek*2
                end
                #
                #to modify Ruser
                Ruser_m[:,dofk]=Ruser_m[:,dofk]+cnstr[i,3]*Ruser_m[:,dofe]
                DOFreg[dofe,1]=0
            end
        end
        #
        #to eliminate columns from Ruser
        k=0
        for i=1:ndof_m
            if DOFreg[i,1]==1
                k=k+1
                Ru[:,k]=Ruser_m[:,i]
            end
        end
        Ruser_m=[]
        Ruser_m=Ru[:,1:k]
        Ruser[(ml-1)*ndof_m+1:ml*ndof_m,(ml-1)*k+1:ml*k]=Ruser_m
    end

    return Ruser

end



function spring_klocal(ku,kv,kw,kq,a,BC,m_a,discrete,ys)
    #
    #Generate spring stiffness matrix [k] in local coordinates, modified from
    #klocal
    #BWS DEC 2015
    
    # Inputs:
    # ku;kv;kw;kq spring stiffness values 
    # a: length of the strip in longitudinal direction
    # BC: ["S-S"] a string specifying boundary conditions to be analyzed:
        #"S-S" simply-pimply supported boundary condition at loaded edges
        #"C-C" clamped-clamped boundary condition at loaded edges
        #"S-C" simply-clamped supported boundary condition at loaded edges
        #"C-F" clamped-free supported boundary condition at loaded edges
        #"C-G" clamped-guided supported boundary condition at loaded edges
    # m_a: longitudinal terms [or half-wave numbers] for this length()
    # discrete .== 1 if discrete spring
    # ys = location of discrete spring
    
    # Output:
    # k: local stiffness matrix; a totalm x totalm matrix of 8 by 8 submatrices.
    # k=[kmp]totalm x totalm block matrix
    # each kmp is the 8 x 8 submatrix in the DOF order [u1 v1 u2 v2 w1 theta1 w2 theta2]'
    
    #
    totalm = length(m_a); #Total number of longitudinal terms m
    #
    #modified by CDM
    k=sparse(zeros(Float64, (8*totalm,8*totalm)))

    z0=zeros(Float64, (4,4))
    for m=1:1:totalm
        for p=1:1:totalm
            #
            km_mp=zeros(Float64, (4,4))
            kf_mp=zeros(Float64, (4,4))
            um=m_a[m]*pi
            up=m_a[p]*pi
            #
            if Int(discrete) == 1
                I1,I5 = BC_I1_5_atpoint(BC,m_a[m],m_a[p],a,ys)
            else #foundation spring
                I1,I2,I3,I4,I5 = BC_I1_5(BC,m_a[m],m_a[p],a)
            end
            #
            #asemble the matrix of km_mp
            km_mp=[ ku*I1   0                   -ku*I1      0
                    0       kv*I5*a^2/(um*up)   0           -kv*I5*a^2/(um*up) 
                    -ku*I1  0                   ku*I1       0
                    0       -kv*I5*a^2/(um*up)  0           kv*I5*a^2/(um*up)]
            #
            #asemble the matrix of kf_mp 
             kf_mp=[ kw*I1   0       -kw*I1      0
                    0       kq*I1   0           -kq*I1 
                    -kw*I1  0       kw*I1       0
                    0       -kq*I1  0           kq*I1];      
            #assemble the membrane & flexural stiffness matrices      
            kmp=[km_mp  z0
                   z0  kf_mp]
    
            #add it into local element stiffness matrix by corresponding to m
            k[8*(m-1)+1:8*m,8*(p-1)+1:8*p] .= kmp
        end
    end
    
    return k    
        
            
    end


function spring_trans(alpha,ks,m_a)
    #
    # Transfer the local stiffness into global stiffness
    # Zhanjie 2008
    # modified by Z. Li; Aug. 09; 2009
    # adapted for spring Dec 2015
    
    totalm = length(m_a) #Total number of longitudinal terms m
    a=alpha
    
    #added by CDM
    gamma = zeros(Float64, (size(ks)[1], size(ks)[2]))
    
    #
    z0 = 0.0
    gam=[cos(a)   z0    z0	  z0	-sin(a)	z0   z0	   z0
            z0	  1	    z0	  z0	  z0    z0	 z0    z0
            z0	  z0  cos(a)  z0	  z0	z0 -sin(a) z0
            z0	  z0	z0	   1	  z0	z0	 z0	   z0
            sin(a) z0	z0	  z0	cos(a)	z0	 z0	   z0
            z0	  z0	z0	  z0	  z0	1	 z0    z0
            z0	  z0  sin(a)  z0	  z0	z0	cos(a) z0
            z0	  z0	z0	  z0	  z0	z0	 z0	   1 ]
    #extend to multi-m
    for i=1:totalm
        gamma[8*(i-1)+1:8*i,8*(i-1)+1:8*i] .= gam
    end
    #
    ksglobal=gamma*ks*gamma'
    
    return ksglobal

end


function spring_assemble(K,k,nodei,nodej,nnodes,m_a)
    #
    #Add the [spring] contribution to the global stiffness matrix

    #Outputs:
    # K: global elastic stiffness matrix
    # K & Kg: totalm x totalm submatrices. Each submatrix is similar to the
    # one used in original CUFSM for single longitudinal term m in the DOF order
    #[u1 v1...un vn w1 01...wn 0n]m'.

    # Z. Li; June 2008
    # modified by Z. Li; Aug. 09; 2009
    # Z. Li; June 2010
    # adapted for springs BWS Dec 2015

    #added by CDM
    nodei = Int(nodei)
    nodej = Int(nodej)
    nnodes = Int(nnodes)

    totalm = length(m_a) #Total number of longitudinal terms m
    K2=sparse(zeros(Float64, (4*nnodes*totalm,4*nnodes*totalm)))
    K3=sparse(zeros(Float64, (4*nnodes*totalm,4*nnodes*totalm)))
    skip=2*nnodes

    for i=1:1:totalm
        for j=1:1:totalm
            #Submatrices for the initial stiffness
            k11=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+1:8*(j-1)+2]
            k12=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+3:8*(j-1)+4]
            k13=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+5:8*(j-1)+6]
            k14=k[8*(i-1)+1:8*(i-1)+2,8*(j-1)+7:8*(j-1)+8]
            k21=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+1:8*(j-1)+2]
            k22=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+3:8*(j-1)+4]
            k23=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+5:8*(j-1)+6]
            k24=k[8*(i-1)+3:8*(i-1)+4,8*(j-1)+7:8*(j-1)+8]
            k31=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+1:8*(j-1)+2]
            k32=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+3:8*(j-1)+4]
            k33=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+5:8*(j-1)+6]
            k34=k[8*(i-1)+5:8*(i-1)+6,8*(j-1)+7:8*(j-1)+8]
            k41=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+1:8*(j-1)+2]
            k42=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+3:8*(j-1)+4]
            k43=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+5:8*(j-1)+6]
            k44=k[8*(i-1)+7:8*(i-1)+8,8*(j-1)+7:8*(j-1)+8]
            #
            K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k11
            if nodej != 0
                K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k12
                K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k21
                K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k22
            end
            #
            K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k33
            if nodej != 0
                K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k34
                K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k43
                K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k44
            end        #
            K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k13
            if nodej != 0
                K2[4*nnodes*(i-1)+nodei*2-1:4*nnodes*(i-1)+nodei*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k14
                K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+skip+nodei*2-1:4*nnodes*(j-1)+skip+nodei*2]=k23
                K2[4*nnodes*(i-1)+nodej*2-1:4*nnodes*(i-1)+nodej*2,4*nnodes*(j-1)+skip+nodej*2-1:4*nnodes*(j-1)+skip+nodej*2]=k24
            end
            #
            K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k31
            if nodej != 0
                K2[4*nnodes*(i-1)+skip+nodei*2-1:4*nnodes*(i-1)+skip+nodei*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k32
                K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+nodei*2-1:4*nnodes*(j-1)+nodei*2]=k41
                K2[4*nnodes*(i-1)+skip+nodej*2-1:4*nnodes*(i-1)+skip+nodej*2,4*nnodes*(j-1)+nodej*2-1:4*nnodes*(j-1)+nodej*2]=k42
            end
            #
        end
    end
    K=K+K2

    return K

end



function strip(prop,node,elem,lengths,springs,constraints,neigs)
    #HISTORY
    #June 2010; complete update to new Boundary conditions; Z. Li; B. Schafer
    #
    #INPUTS
    #prop: [matnum Ex Ey vx vy G] 6 x nmats
    #node: [node# x z dofx dofz dofy dofrot stress] nnodes x 8
    #elem: [elem# nodei nodej t matnum] nelems x 5
    #lengths: [L1 L2 L3...] 1 x nlengths; lengths to be analyzed
    #could be half-wavelengths for signiture curve
    #or physical lengths for general b.c.
    #springs: [node# d.o.f. kspring kflag] where 1=x dir 2= z dir 3 = y dir 4 = q dir (twist) flag says if k is a foundation stiffness | a total stiffness
    #constraints:: [node#e dofe coeff node#k dofk] e=dof to be eliminated k=kept dof dofe_node = coeff*dofk_nodek
    #BC: ["S-S"] a string specifying boundary conditions to be analyzed:
    #"S-S" simply-pimply supported boundary condition at loaded edges
    #"C-C" clamped-clamped boundary condition at loaded edges
    #"S-C" simply-clamped supported boundary condition at loaded edges
    #"C-F" clamped-free supported boundary condition at loaded edges
    #"C-G" clamped-guided supported boundary condition at loaded edges
    #m_all: m_all[length#]=[longitudinal_num# ... longitudinal_num#],longitudinal terms m for all the lengths in cell notation
    # each cell has a vector including the longitudinal terms for this length()
    #neigs - the number of eigenvalues to be determined at length (default=10)

    #OUTPUTS
    #curve: buckling curve [load factor] for each length()
    #curve[l] = [ length mode#1
    #             length mode#2
    #             ...    ...
    #             length mode#]
    #shapes = mode shapes for each length()
    #shapes[l] = mode, mode is a matrix, each column corresponds to a mode.

    #FUNCTIONS CALLED IN THIS ROUTINE
    # \analysis\addspring.m : add springs to K
    # \analysis\assemble.m : asseemble global K; Kg
    # \analysis\elemprop.m : element properties
    # \analysis\kglocal.m : element kg matrix
    # \analysis\klocal.m : element k matrix
    # \analysis\trans.m : trasnform k; kg matrix
    # \analysis\msort.m : clean up 0's; multiple longitudinal terms. | out-of-order terms
    # \analysis\constr_BCFlag.m : determine flags for user constraints & internal [at node] B.transpose(C)s
    # \analysis\cFSM\constr_user.m : user defined contraints in cFSM style

    #CDM Notes
    #April 2021
    #Ported to Julia.
    #Old school CUFSM here, no cFSM and just "S-S" boundary conditions.

    m_all = ones(length(lengths))
    BC = "S-S"

    #MATRIX SIZES
    nnodes = length(node[:,1])
    nelems = length(elem[:,1])
    nlengths = length(lengths)

    #CLEAN UP INPUT
    #clean up 0's; multiple terms. | out-of-order terms in m_all
    # m_all=msort(m_all)

    #DETERMINE FLAGS FOR USER CONSTRAINTS AND INTERNAL [AT NODE] B.transpose(C)s
    BCFlag=constr_BCFlag(node,constraints)

    #GENERATE STRIP WIDTH AND DIRECTION ANGLE
    elprop=elemprop(node,elem,nnodes,nelems)

    #---------------------------------------------------------------------------------

    #added by CDM
    #mimic the Matlab cell structure
    curve = Vector{Matrix{Float64}}(undef, nlengths)
    shapes = Vector{Matrix{Float64}}(undef, nlengths)

    #LOOP OVER ALL THE LENGTHS TO BE INVESTIGATED
    l=0; #length_index = one
    while l<nlengths
        l=l+1; #length index = length index + one

        #length to be analyzed
        a=lengths[l]
        #longitudinal terms to be included for this length()
        m_a=m_all[l]
        #
        totalm=length(m_a);#Total number of longitudinal terms

        #ZERO OUT THE GLOBAL MATRICES
        K=sparse(zeros(Float64, (4*nnodes*totalm,4*nnodes*totalm)))
        Kg=sparse(zeros(Float64, (4*nnodes*totalm,4*nnodes*totalm)))
        #
        #ASSEMBLE THE GLOBAL STIFFNESS MATRICES
        for i=1:nelems
            #Generate element stiffness matrix [k] in local coordinates
            t=elem[i,4]
            b=elprop[i,2]
            matnum=elem[i,5]
            row=findfirst(x->x==matnum, prop[:,1])  #revised by CDM
            Ex=prop[row,2]
            Ey=prop[row,3]
            vx=prop[row,4]
            vy=prop[row,5]
            G=prop[row,6]
            k_l=klocal(Ex,Ey,vx,vy,G,t,a,b,BC,m_a)
            #Generate geometric stiffness matrix [kg] in local coordinates
            Ty1=node[Int(elem[i,2]),8]*t
            Ty2=node[Int(elem[i,3]),8]*t
            kg_l=kglocal(a,b,Ty1,Ty2,BC,m_a)
            #Transform k & kg into global coordinates
            alpha=elprop[i,3]
            k,kg=trans(alpha,k_l,kg_l,m_a)
            #Add element contribution of k to full matrix K & kg to Kg
            nodei=Int(elem[i,2])
            nodej=Int(elem[i,3])
            K,Kg=assemble(K,Kg,k,kg,nodei,nodej,nnodes,m_a)
        end

        #ADD SPRING CONTRIBUTIONS TO STIFFNESS
        #Prior to version 4.3 this was the springs method
            #     if ~isempty(springs) #springs variable exists
            #         [K]=addspring[K,springs,nnodes,a,BC,m_a]
            #     end
        #Now from version 4.3 this is the new springs method
        if (!isempty(springs))
            if (length(springs[1,:])==10) & (springs[1,1]!==0) #springs variable exists, is right length, & non-zero
                nsprings = length(springs[:,1])
                for i=1:nsprings
                    #Generate spring stiffness matrix [ks] in local coordinates
                    ku=springs[i,4]
                    kv=springs[i,5]
                    kw=springs[i,6]
                    kq=springs[i,7]
                    discrete=springs[i,9]
                    ys=springs[i,10]*a
                    ks_l=spring_klocal(ku,kv,kw,kq,a,BC,m_a,discrete,ys)
                    #Transform ks into global coordinates
                    nodei = springs[i,2]
                    nodej = springs[i,3]
                    if nodej==0 #spring is to ground
                        #handle the spring to ground during assembly
                        alpha=0; #use global coordinates for spring
                    else #spring is between nodes
                        xi = node[nodei,2]
                        zi = node[nodei,3]
                        xj = node[nodej,2]
                        zj = node[nodej,3]
                        dx = xj - xi
                        dz = zj - zi
                        width = sqrt(dx^2 + dz^2)
                        if width<1E-10 #coincident nodes
                            alpha=0; #use global coordinates for spring
                        elseif springs[i,8]==0 #flagged for global()
                            alpha=0; #use global coordinates for spring
                        else
                            alpha = atan2(dz,dx); #local orientation for spring
                        end
                    end
                    ks=spring_trans(alpha,ks_l,m_a)
                    #Add element contribution of ks to full matrix K
                    K=spring_assemble(K,ks,nodei,nodej,nnodes,m_a)
                end
            end
        end


    #INTERNAL BOUNDARY CONDITIONS [ON THE NODES] AND USER DEFINED CONSTR.
        #Check for user defined constraints too
        if BCFlag==0
            #no user defined constraints & fixities.
            Ru0=0
            nu0=0
        else
            #size boundary conditions & user constraints for use in R format()
            #d_constrained=Ruser*d_unconstrained; d=nodal DOF vector (note by
            #BWS June 5 2006)
            Ruser=constr_user(node,constraints,m_a)
            Ru0=nullspace(Ruser')
            #Number of boundary conditiions & user defined constraints = nu0
            nu0=length(Ru0[1,:])
        end

        #modified by CDM
        #no modal constraints are activated therefore
        Rmode = Matrix(I, 4*nnodes*totalm, 4*nnodes*totalm)  #activate modal constraints



        #
        #CREATE FINAL CONSTRAINT MATRIX
        #Determine the number of modal constraints; nm0
        if BCFlag==0
            #if no user defined constraints & fixities.
            R=Rmode
        else
            R=nullspace(Ru0')
        end
        #
        #INTRODUDCE CONSTRAINTS AND REDUCE K MATRICES TO FREE PARTS ONLY
        Kff=R'*K*R
        Kgff=R'*Kg*R

        #SOLVE THE EIGENVALUE PROBLEM

        lf=eigvals(Kff,Kgff)
        modes=eigvecs(Kff,Kgff)

        #CLEAN UP THE EIGEN SOLUTION
        #find all the positive eigenvalues & corresponding vectors; squeeze out the rest
        index=findall(x->x>0, lf)
        lf=lf[index]
        modes=modes[:,index]
        #sort from small to large
        lf=sort(lf)
        index = sortperm(lf)
        modes=modes[:,index]
        #only the real part is of interest [eigensolver may give some small nonzero imaginary parts]
        lf=real.(lf)
        modes=real.(modes)
        #
        #truncate down to reasonable number of modes to be kept
        num_pos_modes=length(lf)
        nummodes=minimum([neigs, num_pos_modes])
        lf=lf[1:nummodes]
        modes=modes[:,1:nummodes]
        #
        #FORM THE FULL MODE SHAPE BY BRINGING BACK ELIMINATED DOF
        mode=R*modes
        #
        #CLEAN UP NORMALIZATION OF MODE SHAPE
        #eig & eigs solver use different normalization
        #set max entry [absolute] to +1.0 & scale the rest
        for j=1:nummodes
            maxindex=findall(x->x==maximum(abs.(mode[:,j])), abs.(mode[:,j]))
            mode[:,j]=mode[:,j] ./ mode[maxindex[1],j]
        end
        #

        #added by CDM
        curve[l] = Matrix{Float64}(undef, nummodes, 2)
        shapes[l] = Matrix{Float64}(undef, nnodes*4, nummodes)
        #GENERATE OUTPUT VALUES
        #curve & shapes are changed to cells!!
        #curve: buckling curve [load factor]
        #curve[l] = [ length ... length()
        #             mode#1 ... mode#]
        #shapes = mode shapes
        #shapes[l] = mode, each column corresponds to a mode.
        curve[l][1:nummodes,1] .= lengths[l]
        curve[l][1:nummodes,2] = lf
        #shapes[:,l,1:min([nummodes,num_pos_modes])]=modes
        shapes[l] = mode

    end
    #THAT ENDS THE LOOP OVER ALL THE LENGTHS
    #--------------------------------------------------------------------------

    return curve,shapes
end


function stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)
    #BWS
    #1998
    if unsymm==0
       Ixz=0
    end
    stress=zeros(length(node[:,1]))
    #from P
    stressinc=P/A
    if maximum(isnan.(stressinc))
        stressinc=0.0
    end
    stress = stress .+ stressinc
    #from Mxx Mzz
    stressinc=-((Mzz*Ixx+Mxx*Ixz) .* (node[:,2] .-xcg) .-(Mzz*Ixz+Mxx*Izz) .*(node[:,3] .-zcg)) ./(Izz*Ixx-Ixz^2)
    if maximum(isnan.(stressinc))
        stressinc=0.0
    end
    stress = stress .+ stressinc
    #or from M11 & M22
    th=thetap*pi/180
    prin_coord=inv([cos(th) -sin(th);  sin(th) cos(th)])*[(node[:,2] .-xcg)' ; (node[:,3] .-zcg)']
    #M11
    stressinc = M11 .* prin_coord[2,:] ./I11
    if maximum(isnan.(stressinc))
        stressinc = 0.0
    end
    stress = stress .+ stressinc
    #M22
    stressinc=-M22 .* prin_coord[1,:] ./I22
    if maximum(isnan.(stressinc))
        stressinc = 0.0
    end
    stress = stress .+ stressinc
    #assign stress to node variable
    node[:,8] .= stress

    return node

end



end #module
