"""
    surrogate_model(plan::AbstractArray{T,2}, samples::AbstractArray{T,2}; options=Options()) where T

Returns a surrogate model function based on an optimized Radial Basis Function
interpolant. Depending on the supplied options; the kernel, kernel width and scaling of 
input data is optimized.

...
# Arguments
- `plan::AbstractArray{T,2}`: 
    sample locations where each column corresponds to the location of one point.
    `size(plan) = (num_dimensions,num_samples)`.
- `samples::AbstractArray{T,2}`:
    function value at each sample location. each column contains one value from
    the corresponding plan location. `size(samples) = (1,num_samples)`.
- `options=Options()`: 
    all options available to customize the surrogate optimization.  
...
"""
function surrogate_model(plan::AbstractArray{T,2}, samples::AbstractArray{T,2}; options::Options=Options()) where T

    @unpack num_interpolants, trace, parallel_surrogate, smooth = options

    (trace == :verbose) && println("Creating optimized surrogate model")

    (length(samples) != size(plan,2)) && error("plan and samples do not have the correct sizes")
    
    #Optimize RBF hypers for the ensamble of interpolants
    if parallel_surrogate
        mres = pmap(  (x)->rbf_hypers_opt(samples, plan, options), 
                        1:num_interpolants)
    else
        mres = map(  (x)->rbf_hypers_opt(samples, plan, options), 
                        1:num_interpolants)
    end

    optres = first.(mres)
    fitness = last.(mres)

    #Print hyperparameter results
    (trace == :verbose) && print_sm_opt(fitness,optres,smooth)

    # Sample order for ScatteredInterpolation
    samples_vec = vec(samples)

    # Interpolation object based on the optimization results
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

function print_sm_opt(fitness,optres,smooth)
    print(@sprintf("    mean surrogate fitness (smaller is better)= %.5g,  (min = %.5g, max = %.5g)\n",mean(fitness),minimum(fitness),maximum(fitness)))

    sc_res = [opt.scaling./opt.scaling[1] for opt in optres]
    println("    mean dim. scale (norm. to first dim.) = ", mean(sc_res),
                                        " (min = ",minimum(sc_res),
                                        ", max = ",maximum(sc_res),")")
    sm_res = [opt.smooth for opt in optres]
    if smooth == :variable
        println("    mean smoothing = ", mean.(sm_res),
                                        " (min = ",minimum.(sm_res),
                                        ", max = ",maximum.(sm_res),") \n")
    else
        println(@sprintf("    mean smoothing = %.5g,  (min = %.5g, max = %.5g)\n",mean(sm_res),minimum(sm_res),maximum(sm_res)))
    end
end