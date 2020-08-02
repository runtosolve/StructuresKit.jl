# StructuresKit

[![Build Status](https://ci.appveyor.com/api/projects/status/github/runtosolve/StructuresKit.jl?svg=true)](https://ci.appveyor.com/project/runtosolve/StructuresKit-jl)
[![Codecov](https://codecov.io/gh/runtosolve/StructuresKit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/runtosolve/StructuresKit.jl)

## Usage
This package provides coded pieces and parts for analyzing structural components and systems.

##Modules

[PlautBeam](https://github.com/runtosolve/StructuresKit.jl/blob/master/docs/PlautBeam.md)

Perform second order structural analysis of single or multi-span thin-walled beams with a uniform loading.   The location of the uniform loading on the cross-section (e.g., top of flange, bottom of flange, through the shear center) can be specified.  Lateral and vertical loads can be applied. Continuous lateral and rotational springs are available to simulate attachment bracing.

[InternalForces](https://github.com/runtosolve/StructuresKit.jl/blob/master/docs/InternalForces.md)

Calculate internal axial force, shear, moment, torsion, and bimoment from a structural member's displaced shape.

##Example

```
using StructuresKit



```
