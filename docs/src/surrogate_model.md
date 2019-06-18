# Surrogate model
The surrogate model which is used for the optimisation can be created manually with

```@docs
surrogate_model
```

This enables the freedom to chose how it is optimised and used. When called, the 
function value from each surrogate in the ensemble is returned.

## Example
```julia
julia> rosenbrock_2D(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
julia> search_range=[(-5.0,5.0),(-5.0,5.0)]

# Start from 5 Latin Hypercube Samples
julia> num_samples=5
julia> sampling_plan_opt_gens=10_000
julia> plan = scaled_LHC_sampling_plan(search_range,num_samples,sampling_plan_opt_gens;trace=false)

# Evaluate the function 
julia> samples = mapslices(rosenbrock_2D,plan,dims=1)

# Create the optimized surrogate model (optres contains the optimisation results for the surrogate)
opt=SurrogateModelOptim.Options()
julia> sm_interpolant, optres = surrogate_model(plan, samples;options=opt)
```

# Model infill
New design locations can be found with
```@docs
model_infill
```
 
## Example
```julia
julia> infill_plan, infill_type, infill_prediction, options = model_infill(search_range,plan,samples,sm_interpolant; options=Options())
```
!!! note

    The options are updated to cycle through the `infill_type`.