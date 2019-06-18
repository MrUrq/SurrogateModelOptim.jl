<img src="docs/src/assets/logo.png" width="180">

# SurrogateModelOptim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MrUrq.github.io/SurrogateModelOptim.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://MrUrq.github.io/SurrogateModelOptim.jl/latest)
[![Build Status](https://travis-ci.org/MrUrq/SurrogateModelOptim.jl.svg?branch=master)](https://travis-ci.org/MrUrq/SurrogateModelOptim.jl)
[![Codecov](https://codecov.io/gh/MrUrq/SurrogateModelOptim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MrUrq/SurrogateModelOptim.jl)
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- [![Build Status](https://ci.appveyor.com/api/projects/status/github/MrUrq/SurrogateModelOptim.jl?svg=true)](https://ci.appveyor.com/project/MrUrq/SurrogateModelOptim-jl) -->

<!-- [![Coveralls](https://coveralls.io/repos/github/MrUrq/SurrogateModelOptim.jl/badge.svg?branch=master)](https://coveralls.io/github/MrUrq/SurrogateModelOptim.jl?branch=master) -->


*SurrogateModelOptim* is a Julia package for the optimisation of expensive functions. 
The surrogate model is based on an ensemble of Radial Basis Function interpolants with adaptive axis scaling.

Features:

* Sampling plan creation through 
* Creation of an optimised RBF surrogate.
* Infill of design space.

## Installation

The package is registered and can be installed with `Pkg.add`.

```julia
julia> Pkg.add("SurrogateModelOptim")
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **tagged version of the documentation.**


## Author

- Magnus Urquhart - [@MrUrq](https://github.com/MrUrq/)

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MrUrq.github.io/SurrogateModelOptim.jl/stable

<!-- ### Reference
This package is a 
[1]: Stuart Bates, Johann Sienz, and Vassili Toropov. "Formulation of the Optimal Latin Hypercube Design of Experiments Using a Permutation Genetic Algorithm", 45th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics & Materials Conference, Structures, Structural Dynamics, and Materials and Co-located Conferences, () https://doi.org/10.2514/6.2004-2011 -->