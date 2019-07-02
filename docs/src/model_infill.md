# Model infill
New design locations can be found with
```@docs
model_infill
```
 
The infill criteria is optimized based on the supplied options. Custom criteria can be
implemented and solved for by using the `surrogate_model` with your favourite optimization
library.  

## Example
The package supplied model infill can be used as 
```julia
julia> infill_plan, infill_type, infill_prediction, options = model_infill(search_range,plan,samples,sm_interpolant; options=Options())
```
!!! note

    The options are updated to cycle through the `infill_type`.