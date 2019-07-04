# Model infill
New design locations can be found with
```@docs
model_infill
```
 
The infill criteria is optimized based on the supplied options. Custom criteria can be
implemented and solved for by using the `surrogate_model` with your favourite optimization
library.  

In the case where several infill points are requested, each infill criteria is used
cyclically. If the number of requested infill points exceeds the number of infill
criteria, a point from the pareto front of the infill objectives is selected. The selected
points is the point which maximizes the distance to the already sampled points to explore
the search space. In the case where the requested infill points can't be found without
duplicates, a random point is added. This happens very rarely. Duplicates are never added.
The infill type can be seen during run time with `trace=:verbose`, alternatively the
returned results contain all the infill criteria used. 

## Example
The package supplied model infill can be used as 
```julia
julia> infill_plan, infill_type, infill_prediction, options = model_infill(search_range,plan,samples,sm_interpolant; options=SurrogateModelOptim.Options())
```
!!! note

    The options are updated to cycle through the `infill_type`.