# Surrogate model optimization
A julia function can be optimised with

```@docs
smoptimize
```
 
## Example
```julia
julia> using SurrogateModelOptim
julia> rosenbrock_2D(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
julia> search_range=[(-5.0,5.0),(-5.0,5.0)]
julia> smoptimize(rosenbrock_2D, search_range)
```

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
The default option `num_interpolants=20` meaning the surrogate model ensemble
contains 20 RBF interpolants has shown good performance for a variety of functions.
More test functions and examples can be found in the source files under `examples`.
The performance is typically good for smooth functions with or without noise. 
Discontinuous functions are not captured well by RBF interpolation, see 
the example folder to see this problem.