<img src="docs/src/assets/logo.png" width="180">

# SurrogateModelOptim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MrUrq.github.io/SurrogateModelOptim.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://MrUrq.github.io/SurrogateModelOptim.jl/latest)
[![Build Status](https://travis-ci.org/MrUrq/SurrogateModelOptim.jl.svg?branch=master)](https://travis-ci.org/MrUrq/SurrogateModelOptim.jl)
[![Codecov](https://codecov.io/gh/MrUrq/SurrogateModelOptim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MrUrq/SurrogateModelOptim.jl)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- [![Build Status](https://ci.appveyor.com/api/projects/status/github/MrUrq/SurrogateModelOptim.jl?svg=true)](https://ci.appveyor.com/project/MrUrq/SurrogateModelOptim-jl) -->

<!-- [![Coveralls](https://coveralls.io/repos/github/MrUrq/SurrogateModelOptim.jl/badge.svg?branch=master)](https://coveralls.io/github/MrUrq/SurrogateModelOptim.jl?branch=master) -->


*SurrogateModelOptim* is a Julia package for the optimisation of expensive functions. 
The surrogate model is based on an ensemble of Radial Basis Function interpolants with adaptive axis scaling.

## Installation

The package is registered and can be installed with `Pkg.add`.

```julia
julia> Pkg.add("SurrogateModelOptim")
```

## Optimization
This package is intended to be used for functions which are expensive. Expensive
is in this case considered a function that evaluates in several minutes to days.
The simplest form of usage is as follows.
```julia
julia> using SurrogateModelOptim
julia> rosenbrock_2D(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
julia> search_range=[(-5.0,5.0),(-5.0,5.0)]
julia> smoptimize(rosenbrock_2D, search_range)
```
There are many options accessible through the options interface. The target is to minimize
the function value. The model is created from a Latin Hypercube sampling plan. Several
Radial Basis Function surrogate models are fitted to the data where the ensemble of
surrogates is used to predict new design locations. New designs are added which exploits
the surrogate model as well as explores the design space.

Due to the high cost of creating several surrogates it is highly advisable to create
the surrogate model in parallel. Start julia in parallel with `> julia -p x` where `x`
is the number of available cores. The previous example can then be run as
```julia
julia> result = smoptimize(rosenbrock_2D, search_range;
                    options=SurrogateModelOptim.Options(
                    iterations=25,
                    num_interpolants=N*x, #Where N is an integer number
                    num_start_samples=5,
                        ));
```
`num_interpolants=10` meaning that the surrogate model ensemble contains 10 RBF interpolants
 has shown good performance for a variety of functions. 

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **tagged version of the documentation.**


## Author

- Magnus Urquhart - [@MrUrq](https://github.com/MrUrq/)

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MrUrq.github.io/SurrogateModelOptim.jl/stable

### Citation
TBD
