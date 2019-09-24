# Examples
Most of the flexibility is left for the user to create and tailor to their specific needs.
There are several examples located in `examples/` of different ways the optimization
package can be used. To run all simply do
```julia
julia> include(joinpath(dirname(pathof(SurrogateModelOptim)), "../examples/test_all.jl"))
```
`BlackBoxOptim.jl, PlotlyJS.jl, Distances.jl, LatinHypercubeSampling.jl and Parameters.jl` are needed to run
all examples. Below is a brief explanation of the examples. 

## Noisy optimization
Optimization of noisy functions requires no special treatment except that `smooth` is used
in the options. For most functions, even those with small amounts of noise, it is
recommended to use `smooth=:single`.

## Categorical optimization
The example shows how it is possible to constrain the optimization to discrete values.

## Multi-objective optimization
Multi-objective optimization can be performed by creating a separate surrogate of each
objective. The surrogate can then be used with you favourite algorithm, in this example
[BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) is used. In this
example a random point is selected from the pareto front for infill. It would likely be
better to add a space filling criteria when adding the infill points.

## Discontinuous optimization
This example shows how the performance of the optimization suffers in regions of
discontinuity. If possible, reformulate the problem so that it becomes smooth. If that is
not possible it is better to use another optimization package, especially if the optimum
is expected close to the discontinuity. Note that noise can be handled and is not
considered to be a discontinuity. 

## Constrained optimization
Constrained optimization can be performed by adding an additional penalty after the
surrogate model has been created. Note that it is important to add it after the creation
to keep the function as smooth as possible. If it is added before the creation, the
surrogate accuracy will suffer in the transition area where the penalty was added.

## Fast options
Run an example using the options for running the optimisation fast. This is not
a recommended settings to use. Use another optimisation package such as `Optim.jl` or
`BlackBoxOptim.jl` if the function evaluation is fast compared to the surrogate model creation.
Can be useful for faster validation of code. 