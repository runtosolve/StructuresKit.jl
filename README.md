# StructuresKit.jl

*Analyze and design structural systems* 


| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |


## Installation

This package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add StructuresKit
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("StructuresKit")
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*

## Project Status

This package is under heavy development.   The plan is to offer a set of primitives (building blocks) like `beam, column, shell, connection` and then make it fast and easy to assemble, analyze, design, and visualize the structural system.   

## Questions and Contributions

Usage questions can be posted in the #structures-kit channel on [Julia Slack](https://julialang.org/community/).

Contributions are very welcome, as are feature requests and suggestions. Please open an [issue][issues-url] if you encounter problems. Please open pull requests against the `master` branch whenever possible. 


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://github.com/runtosolve/StructuresKit.jl/tree/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://github.com/runtosolve/StructuresKit.jl/tree/master

[travis-img]: https://travis-ci.org/runtosolve/StructuresKit.jl.svg?branch=master
[travis-url]: https://travis-ci.org/StructuresKit.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/idfm6woehn70umgn?svg=true
[appveyor-url]: https://ci.appveyor.com/project/cristophermoen/structureskit-jl

[codecov-img]: https://codecov.io/gh/runtosolve/StructuresKit.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/runtosolve/StructuresKit.jl

[issues-url]: https://github.com/runtosolve/StructuresKit.jl/issues


