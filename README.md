# StructuresKit

[![Build Status](https://travis-ci.org/runtosolve/StructuresKit.jl.svg?branch=master)](https://travis-ci.org/runtosolve/StructuresKit.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/idfm6woehn70umgn?svg=true)](https://ci.appveyor.com/project/cristophermoen/structureskit-jl)
[![codecov](https://codecov.io/gh/runtosolve/StructuresKit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/runtosolve/StructuresKit.jl)


## Usage
Analyze and design structural components and systems with a variety of Julia modules.   

## Status
The current offerings are pretty niche, with a focus on thin-walled structures.  There is plenty of potential for growth, especially for large system problems. This is just the beginning.

## Modules

### Analysis

[PlautBeam](https://github.com/runtosolve/StructuresKit.jl/blob/master/docs/PlautBeam/PlautBeam.md)

Perform second order structural analysis of single or multi-span thin-walled beams with a uniform loading.  

[InternalForces](https://github.com/runtosolve/StructuresKit.jl/blob/master/docs/InternalForces/InternalForces.md)

Calculate internal axial force, shear, moment, and torsion from a structural member's displaced shape.

### Codes and Standards

[AISIS10016](https://github.com/runtosolve/StructuresKit.jl/blob/master/docs/AISIS10016/AISIS10016.md)

American Iron and Steel Institute (AISI) S100-16 *North American Specification for the Design of Cold-Formed Steel Structural Members*

[AISIS10024](https://github.com/runtosolve/StructuresKit.jl/blob/master/docs/AISIS10024/AISIS10024.md)

American Iron and Steel Institute (AISI) S100-24 *North American Specification for the Design of Cold-Formed Steel Structural Members*

## Version planning

0.1  current
0.2  add PlautColumn
0.X  add LocalBuckling
0.x  add Dynamics
0.x  add Visualization
1.0  upgrade to matrix free methods
2.0  add member connectivity
3.0  add multi-threading and parallel
