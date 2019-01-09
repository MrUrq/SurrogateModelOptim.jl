"""
    function interp_obj(inpt::Vector{Float64}, kerns, samples::Array{Float64,2},
    plan::Array{Float64,2}; rippa::Bool = false, scale::Bool = false)

Objective function for optimisation of interpolation kernel function and width.
"""
function interp_obj(inpt::Vector{Float64}, kerns, samples,
        plan::Array{Float64,2}; rippa::Bool = false,
        variable_kernel_width::Bool = true, variable_dim_scaling::Bool = true,
        smooth = false, cond_max=cond_max, rbf_dist_metric = Distances.Euclidean(),
        smooth_user::Float64 = 0.0)

    
    kern, scaling, smooth = extract_bboptim_hypers( inpt,plan,kerns,variable_kernel_width,
                                                    variable_dim_scaling,smooth,
                                                    smooth_user)
    
    optres = RBFHypers(kern, scaling, smooth)

    #Preprocess the plan based on the settings used.
    preprocessed_plan = preprocess_point(plan,optres,base_scale=plan)

    
    E = try
        E = RMSErrorLOO(kern, samples, preprocessed_plan, smooth; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
    catch 
        E = Inf
    end

    return E
end


function rbf_hypers_opt(samples_org::Array{Float64,2}, plan::Array{Float64,2}, options::SurrogateModelOptim.Options)
    
    @unpack rippa, variable_kernel_width, variable_dim_scaling, rbf_opt_method, 
            min_rbf_width, max_rbf_width, min_scale, max_scale, cond_max,
            rbf_dist_metric, rbf_opt_gens, kerns, smooth, max_smooth, smooth_user = options

    samples = vec(samples_org)


    # Create the hyperparameter search range based on the input options 
    sr = construct_search_range(plan, variable_kernel_width,
                                min_rbf_width, max_rbf_width, variable_dim_scaling,
                                min_scale, max_scale, smooth, max_smooth)

    # RBF hyperparameter objective function
    itp_obj = function (x)
        interp_obj(x,kerns,samples,plan; 
                rippa=rippa, variable_kernel_width=variable_kernel_width,
                variable_dim_scaling=variable_dim_scaling, smooth=smooth,
                cond_max=cond_max)
    end

    # Optimize the interpolant hyperparameters
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:silent, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6);
        
    kern, scaling, smooth = extract_bboptim_hypers( res.archive_output.best_candidate,
                                                    plan,kerns,variable_kernel_width,
                                                    variable_dim_scaling,smooth,smooth_user)

    # Return the optimized hyperparameters in the correct type
    return RBFHypers(kern, scaling, smooth)
end





#Returns anonymous function that is an optimised RBF-based surrogate model
function surrogate_model(samples, plan, options)
    
    @unpack num_interpolants = options
    
    #Optimize RBF hypers for the ensamble of interpolants
    optres = pmap(  (x)->rbf_hypers_opt(samples, plan, options), 
                    1:num_interpolants)

    #Correct sample order
    observations = samples'

    #Interpolation object based on the optimisation results
    itp = Array{ScatteredInterpolation.RBFInterpolant,1}(undef,num_interpolants)
    for i = 1:num_interpolants
        #Preprocess the points based on the settings used.
        preprocessed_point = preprocess_point(plan,optres[i],base_scale=plan)
        
        itp[i] = interpolate(optres[i].kernelFunc, preprocessed_point,
                             observations, smooth=optres[i].smooth)
    end

    return function (estimation_point)        
        res = Array{Float64,2}(undef,num_interpolants,1)

        for i = 1:num_interpolants
            preprocessed_est_point = preprocess_point(estimation_point,optres[i],base_scale=plan)

            #Evaluate the interpolation object
            res[i] = ScatteredInterpolation.evaluate(itp[i], preprocessed_est_point)[1]
        end
        return SurrogateEstimate(res)        
    end
end



