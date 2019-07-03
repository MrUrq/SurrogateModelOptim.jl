# Surrogate model
The surrogate model which is used for the optimization can be created manually with

```@docs
surrogate_model
```

This enables the freedom to choose how it is optimized and used. When called, the 
function value from each surrogate in the ensemble is returned.

## Example
```julia
julia> rosenbrock_2D(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
julia> search_range=[(-5.0,5.0),(-5.0,5.0)]

# Start from 5 Latin Hypercube Samples
julia> num_samples=5
julia> sampling_plan_opt_gens=10_000
julia> plan = scaled_LHC_sampling_plan(search_range,num_samples,sampling_plan_opt_gens)

# Evaluate the function 
julia> samples = mapslices(rosenbrock_2D,plan,dims=1)

# Create the optimized surrogate model (optres contains the optimization results for the surrogate)
julia> opt=SurrogateModelOptim.Options()
julia> sm_interpolant, optres = surrogate_model(plan, samples;options=opt)
```

The surrogate model is by default created using only Gaussian kernels. This can be changed
by supplying several kernels, e.g. `kerns =
[ScatteredInterpolation.Gaussian,
ScatteredInterpolation.InverseQuadratic,
ScatteredInterpolation.InverseMultiquadratic]`.
For a limited number of `rbf_opt_gens` the surrogate model performance tends to be better
using only one kernel to limit the search space when optimizing the surrogate
hyperparameters. 

Turning on verbosity with `trace=:verbose` can help with judging of the hyperparameter
optimisation.