"""
    surrogate_model(plan::AbstractArray{T,2}, samples::AbstractArray{T,2}; options=Options()) where T

Returns surrogate model that is a function based on an optimised Radial Basis Function
interpolant. Depending on the options, the kernel, kernel width and scaling of 
input data is optimised.

...
# Arguments
- `plan::AbstractArray{T,2}`: 
    sample locations where each column corresponds to the location of one point.
    `size(plan) = (num_dimensions,num_samples)`.
- `samples::AbstractArray{T,2}`:
    function value at each sample location. each column contains one value from
    the corresponding plan location. `size(samples) = (1,num_samples)`.
- `options=Options()`: 
    all options available to customize the surrogate optimisation.  
...
"""
function surrogate_model(plan::AbstractArray{T,2}, samples::AbstractArray{T,2}; options::Options=Options()) where T

    @unpack num_interpolants, trace, parallel_surrogate = options

    trace && println("Creating optimized surrogate model ...")

    (length(samples) != size(plan,2)) && error("plan and samples do not have the correct sizes")
    
    #Optimize RBF hypers for the ensamble of interpolants
    if parallel_surrogate
        optres = pmap(  (x)->rbf_hypers_opt(samples, plan, options), 
                        1:num_interpolants)
    else
        optres = map(  (x)->rbf_hypers_opt(samples, plan, options), 
                        1:num_interpolants)
    end

    # Sample order for ScatteredInterpolation
    samples_vec = vec(samples)

    # Interpolation object based on the optimisation results
    itp = Array{ScatteredInterpolation.RBFInterpolant,1}(undef,num_interpolants)
    for i = 1:num_interpolants
        # Preprocess the points based on the settings used.
        preprocessed_point = preprocess_point(plan,optres[i],base_scale=plan)
        
        itp[i] = interpolate(optres[i].kernelFunc, preprocessed_point,
                            samples_vec, smooth=optres[i].smooth)
    end

    # Pre allocate and calculate to reduce cost of evaluating surrogate model.
    preprocessed_est_point = Array{Float64,2}(undef,size(plan,1),1)
    old_min = minimum(plan,dims=2)
    old_max = maximum(plan,dims=2)

    sm_func = (estimation_point) ->  surrogate_evaluate.(Ref(preprocessed_est_point),
                                        Ref(estimation_point),Tuple(itp),Tuple(optres),
                                        Ref(old_min),Ref(old_max))

    return sm_func, optres
end