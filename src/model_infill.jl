"""
    model_infill(search_range::Vector{Tuple{Float64,Float64}},plan::AbstractArray{T,2},
    samples::AbstractArray{T,2},sm_interpolant; options::Options=Options()) where T

Infill function that calculates the location of new samples based on the supplied options.
The returned options are updated to facilitate cycling through the infill objective
functions.

...
# Arguments
- `search_range::Vector{Tuple{Float64,Float64}`:
    a vector of tuples containing the lower and upper limits
    for each dimension. `length(search_range)` is equal to number
    of dimensions.
- `plan::AbstractArray{T,2}`: 
    sample locations where each column corresponds to the location of one point.
    `size(plan) = (num_dimensions,num_samples)`.
- `samples::AbstractArray{T,2}`:
    function value at each sample location. each column contains one value from
    the corresponding plan location. `size(samples) = (1,num_samples)`.
- `options=Options()`: 
    all options available to customize the surrogate infill.  
...
"""
function model_infill(search_range::Vector{Tuple{Float64,Float64}},plan::AbstractArray{T,2},
        samples::AbstractArray{T,2},sm_interpolant; options::Options=Options()) where T

    @unpack rbf_opt_gens, num_infill_points, trace,
            infill_funcs, infill_iterations = options
        
    (trace == :verbose) && print_infill_head() 

    (length(samples) != size(plan,2)) && error("plan and samples do not have the correct sizes")
    ((length(sm_interpolant(plan[:,1])) == 1) && (any(infill_funcs.==:std))) && error(":std infill criteria can not be used with one interpolant")
    
    # Extract the requested infill objective functions
    infill_obj_fun = infill_objective(sm_interpolant,plan,samples,infill_funcs)

    # Find the infill optimization points
    infill_plan,infill_type,infill_prediction,res_bboptim,options = 
        infill_opt(search_range,infill_iterations,num_infill_points,infill_obj_fun,infill_funcs,plan,sm_interpolant,options)
        
    # Add additional points to fill up the desired number of infill points if required
    infill_plan,infill_type,infill_prediction = 
        infill_add(sm_interpolant,samples,plan,infill_prediction,search_range,infill_plan,
        infill_type,num_infill_points,infill_obj_fun,infill_funcs,res_bboptim)

    (trace == :verbose) && print_infill_tail(infill_type,num_infill_points,infill_plan,infill_prediction)

    return infill_plan, infill_type, infill_prediction, options
end