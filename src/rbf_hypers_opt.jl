"""
    function interp_obj(inpt::Vector{Float64}, kerns, samples::Array{Float64,2},
    plan::Array{Float64,2}; rippa::Bool = false, scale::Bool = false)

Objective function for optimisation of interpolation kernel function and width.
"""
function interp_obj(inpt::Vector{Float64}, kerns, samples,
        plan::Array{Float64,2}; rippa::Bool = false,
        variable_kernel_width::Bool = true, variable_dim_scaling::Bool = true,
        cond_max=cond_max, rbf_dist_metric = Distances.Euclidean())

    if variable_kernel_width
        if variable_dim_scaling
            return _scaled_interp_varK_varW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        else
            return _interp_varK_varW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        end
    elseif !variable_kernel_width
        if variable_dim_scaling
            return _scaled_interp_fixK_fixW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        else
            return _interp_fixK_fixW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        end        
    end
end




function _rbf_hypers_opt(samples_org::Array{Float64,2}, plan::Array{Float64,2}, options::SurrogateModelOptim.Options)
    
    @unpack rippa, variable_kernel_width, variable_dim_scaling, rbf_opt_method, 
            min_rbf_width, max_rbf_width, min_scale, max_scale, cond_max,
            rbf_dist_metric, rbf_opt_gens, kerns, smooth, max_smooth = options

    samples = vec(samples_org)


    # Create the hyperparameter search range based on the input options 
    sr = construct_search_range(plan, variable_kernel_width,
                                min_rbf_width, max_rbf_width, variable_dim_scaling,
                                min_scale, max_scale, smooth, max_smooth)

    # RBF hyperparameter objective function
    itp_obj = function (x)
        interp_obj(x,kerns,samples,plan; 
                rippa=rippa, variable_kernel_width=variable_kernel_width,
                variable_dim_scaling=variable_dim_scaling, cond_max=cond_max)
    end

    # Optimize the interpolant hyperparameters
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:silent, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6);
        

    # Extract and order the results in an Array of RBFHypers
    kern, scaling, smooth = extract_bboptim_hypers( res.archive_output.best_candidate,
                                                    plan,kerns,variable_kernel_width,
                                                    variable_dim_scaling,smooth)

    # Return the optimized hyperparameters in the correct type
    return RBFHypers(kern, scaling, smooth)
end





############################################################################################Move the interpolation out from this function
function _surrogate_interpolant(optres, points, observations, estimation_point)
    
    #Preprocess the points based on the settings used.
    preprocessed_point = preprocess_point(points,optres,base_scale=points)
    preprocessed_est_point = preprocess_point(estimation_point,optres,base_scale=points)


    #Interpolation object based on the optimisation results
    itp = interpolate(optres.kernelFunc, preprocessed_point, observations)
    
    #Evaluate the interpolation object
    estimation = ScatteredInterpolation.evaluate(itp, preprocessed_est_point)[1]

    return estimation
end



#Returns anonymous function 
function surrogate_model(samples, plan, options)
    
    @unpack num_interpolants = options

    _p_rbf_opt = (x)->_rbf_hypers_opt(samples, plan, options)
    
    optres = pmap(_p_rbf_opt,1:num_interpolants)
    
    return function (estimation_point)        
        res = Array{Float64,2}(undef,num_interpolants,1)

        for i = 1:length(optres)
            res[i]= _surrogate_interpolant(
                                            optres[i],
                                            plan,
                                            samples',
                                            estimation_point,
                                            )
        end
        return SurrogateEstimate(res)        
    end
end





