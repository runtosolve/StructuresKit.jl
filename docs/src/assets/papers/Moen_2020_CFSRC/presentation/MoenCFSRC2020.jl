### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ cb011f50-1189-11eb-0dd7-ed7def6268ea
begin
	using Images
	CFSRC = load("CFSRC.png")
end

# ╔═╡ c62118b2-1302-11eb-2835-8dc23025eb2f
using StructuresKit

# ╔═╡ ea974714-1188-11eb-2eab-9598023d0820
url = "https://cfsrc.org/wp-content/uploads/2020/04/cfsrc_short_colloquium_960x540-900x400.png";

# ╔═╡ a7b90454-1189-11eb-0a76-5bf07b38e64d
download(url,"CFSRC.png");

# ╔═╡ b75fcd8c-13b1-11eb-0e6b-614f18f70a73
md""" # Metal building purlin line strength by computation

 Cristopher D. Moen, Ph.D., P.E.


"""

# ╔═╡ 648eface-13b1-11eb-00a4-a3800f8028ef
load("/Users/crismoen/Downloads/RunToSolveLogo_small.png")

# ╔═╡ 69b9fbc2-12f1-11eb-2579-fdb2399c1487
load("/Users/crismoen/Downloads/Butler-DaytonFreight_1360x640.jpg")

# ╔═╡ cb3b85d4-12ea-11eb-023c-558124a94783
md""" ### Calculating the capacity of a metal building roof purlin line is challenging .

- Zee section or Cee section
- screw fastened, standing seam provides partial bracing to the purlin line
- load eccentric from shear center
- roof slope
- insulation
- lap splices
- simple span or continuous span
- uplift or gravity loads
- intermediate bridging, torsion bracing"""

# ╔═╡ b653844e-12f1-11eb-0c09-bb7589224593
md""" ### Thanks to 60+ years of research we have the answers!

- code-base pathways like AISI S100-16 I6.1, similar in Eurocode
- accessible structural analysis software that can handle weird purlins, e.g., new MASTAN2 
- system-level design enabled with Direct Strength Method
- system tests and simulation
- local connection tests, i.e. F-test, rotational restraint,...
"""

# ╔═╡ 36c69d6c-12fe-11eb-1169-d334f22c7e8d
md""" ### The computation workflow is:

- model the purlin line with bracing and bridging and load eccentricities and roof slope and ...
- apply load to the purlin line
- calculate purlin deformations and force, moment, shear, and torsion demands 
- calculate cross-sectional buckling deformation of the free flange 
- compare demands to strength limit state capacities 
- available purlin line strength is the load where a limit state is reached

"""

# ╔═╡ 3ae9166e-12fd-11eb-17c5-65ea1ea6a1f0
load("/Users/crismoen/Downloads/Moen_ID37.png")

# ╔═╡ e2cf0a8c-13cf-11eb-1c03-f77ed9c340b5
md""" ### This workflow is implemented in the open-source software package [StructuresKit.jl](https://github.com/runtosolve/StructuresKit.jl). """

# ╔═╡ 48fc80c2-13d1-11eb-3e19-fb9e4a062801
load("/Users/crismoen/Downloads/StructuresKit.png")

# ╔═╡ 4128b952-13d0-11eb-2e23-ebab51a78b0a
md""" ### Now let's practice. """

# ╔═╡ 241e4886-1303-11eb-16d8-0f4e0e2f6ac7
md""" Compute the purlin line capacity of an actual experiment from [Gao and Moen (2013)]  (https://www.sciencedirect.com/science/article/pii/S0143974X12002489).  Select specimen Z200B-TH25-2, simple span, uplift loading, screw-fastened panel with 25 mm of Thermax rigid board insulation."""

# ╔═╡ b868d0fc-1307-11eb-0439-713d8e56b7ed
load("/Users/crismoen/Downloads/JCSR.png")

# ╔═╡ 2986dfdc-1309-11eb-034c-23e2bf031595
md"""[YouTube]  (https://www.youtube.com/watch?v=_bZg8S_uPLk)"""

# ╔═╡ f1261c3c-1308-11eb-387c-6f85ffbb3cfd
load("/Users/crismoen/Downloads/ztest.png")

# ╔═╡ 3fdc24da-13ad-11eb-2be2-0328cdd1f1f5
load("/Users/crismoen/Downloads/CZ.png")

# ╔═╡ bd7f0896-13ce-11eb-2ebe-ef72663d67b2
md""" Use a 'kitchen sink' interaction equation that is more or less a yield criterion for the free flange-web intersection."""

# ╔═╡ 4d626470-13ad-11eb-0026-c75eddbfd5c2
md""" $\overline{M_x}/M_{ax{\ell}o} + \overline{M_y}/M_{ay{\ell}o} + \overline{M_yf}/M_{ayf{\ell}o} + \overline{B}/B_a \leq 1.15$ """

# ╔═╡ 97293776-130a-11eb-3d5b-37a31c01694b
begin
	ASDorLRFD = 2;   #nominal, no factors
	GravityOrUplift = 1;
	MemberDefinitions = [(7468.0, 74.68, 1, 1, 1, 1, 1)]; 
	SectionProperties = [(6.023814161033781e6, 1.2383858093197288e6, -2.0172242627006555e6, 2039.513635254465, 8.855200586032066e9, 7572.394863421424, 2.91e7, 9.9999999999e10, 9.9999999999e10)];
	CrossSectionDimensions =  [(2.54, 204, 68.9, 26.6, 46, 1, 185.92000000000002, 32)];
	MaterialProperties = [ (203395, 0.3, 420) ];
	LoadLocation = [ (31.179259190771152, 104.71118275417487)];
	BracingProperties = [(6.7, 928.8824383, 7468.0, 7468.0)];
	RoofSlope = 0.0;
	PurlinSpacing = 2070.0 
	EndBoundaryConditions = [1 1];
	Supports = [0 7468.0];
	Bridging = [0 7468.0];
	
	q_test, z, properties, deformation, strengths, forces, interactions, demand_to_capacity = PurlinDesigner.lineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports, Bridging);
end;
	

# ╔═╡ 0f9cf924-130e-11eb-2a59-89eb74c97365
begin
	using Plots
	plot(z, forces.T, legend= :false)
	ylabel!("T [N-mm]")
	xlabel!("z [mm]")
end

# ╔═╡ 72f64fc6-130d-11eb-2385-6db619596723
ComputedPressure = q_test/PurlinSpacing*10^6

# ╔═╡ c47b15c0-130d-11eb-0c97-b31794c867cc
TestPressure = -650.0

# ╔═╡ 68cf72d2-130f-11eb-251c-dda03a98aa6c
begin	
	plot(z, deformation.ϕ, legend= :false)
	ylabel!("v [mm]")
	xlabel!("z [mm]")
end

# ╔═╡ 5dbeec6e-13d3-11eb-227a-416b34cca9df
begin	
	plot(z, deformation.uf, legend= :false)
	ylabel!("u, free flange [mm]")
	xlabel!("z [mm]")
end

# ╔═╡ 8235b9be-13d1-11eb-3ca0-af1d5396aa84
begin	
	plot(z, demand_to_capacity.BT, legend= :false)
	ylabel!("D/C")
	xlabel!("z [mm]")
end

# ╔═╡ 092bfc64-1310-11eb-2595-39db4563b6c3
md""" ### Conclusions

- Structural analysis of metal building purlin lines is challenging.
- Research investment is paying off.
- Purlin line strength can be computed.

"""

# ╔═╡ 6bb5cfc4-1311-11eb-2dd8-2522669269a6
md""" ### Open-Source Software Pep Talk

Open means open!   Move your stuff to Github, learn open-source languages, code more, and contribute to your favorite packages!

"""

# ╔═╡ Cell order:
# ╟─ea974714-1188-11eb-2eab-9598023d0820
# ╟─a7b90454-1189-11eb-0a76-5bf07b38e64d
# ╟─cb011f50-1189-11eb-0dd7-ed7def6268ea
# ╠═b75fcd8c-13b1-11eb-0e6b-614f18f70a73
# ╟─648eface-13b1-11eb-00a4-a3800f8028ef
# ╟─69b9fbc2-12f1-11eb-2579-fdb2399c1487
# ╟─cb3b85d4-12ea-11eb-023c-558124a94783
# ╟─b653844e-12f1-11eb-0c09-bb7589224593
# ╟─36c69d6c-12fe-11eb-1169-d334f22c7e8d
# ╟─3ae9166e-12fd-11eb-17c5-65ea1ea6a1f0
# ╟─e2cf0a8c-13cf-11eb-1c03-f77ed9c340b5
# ╟─48fc80c2-13d1-11eb-3e19-fb9e4a062801
# ╟─4128b952-13d0-11eb-2e23-ebab51a78b0a
# ╟─241e4886-1303-11eb-16d8-0f4e0e2f6ac7
# ╠═b868d0fc-1307-11eb-0439-713d8e56b7ed
# ╟─2986dfdc-1309-11eb-034c-23e2bf031595
# ╟─f1261c3c-1308-11eb-387c-6f85ffbb3cfd
# ╟─3fdc24da-13ad-11eb-2be2-0328cdd1f1f5
# ╟─bd7f0896-13ce-11eb-2ebe-ef72663d67b2
# ╟─4d626470-13ad-11eb-0026-c75eddbfd5c2
# ╠═c62118b2-1302-11eb-2835-8dc23025eb2f
# ╠═97293776-130a-11eb-3d5b-37a31c01694b
# ╠═72f64fc6-130d-11eb-2385-6db619596723
# ╠═c47b15c0-130d-11eb-0c97-b31794c867cc
# ╠═0f9cf924-130e-11eb-2a59-89eb74c97365
# ╠═68cf72d2-130f-11eb-251c-dda03a98aa6c
# ╠═5dbeec6e-13d3-11eb-227a-416b34cca9df
# ╠═8235b9be-13d1-11eb-3ca0-af1d5396aa84
# ╟─092bfc64-1310-11eb-2595-39db4563b6c3
# ╟─6bb5cfc4-1311-11eb-2dd8-2522669269a6
