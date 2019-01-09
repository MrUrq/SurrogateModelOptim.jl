

# function for creating optimised surrogate
function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}, options=options())

    options = _update_options(search_range;options...)

    #Load the optional argument values
    @unpack num_start_samples, show_trace, sampling_plan_opt_gens = options
    
    #Create sampling plan
    plan = _LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens)

    println("Iteration ", 1)
    println("Evaluating function")
    #Evaluate the expensive function
    samples = mapslices(f,plan,dims=1)
    @show minimum(samples)

    sm_interpolant = surrogate_model(samples, plan, options)
        
    ##################################### Bad way of doing this, just proof of concept
    for i = 2:15

        println("Creating surrogate model")
        #Create the optimized Radial Basis Function interpolant 
        sm_interpolant = surrogate_model(samples, plan, options)
        
        println("Finding infill points")
        #Points to add to the sampling plan to improve the interpolant
        plan = model_infill(plan,samples,sm_interpolant,options)
        plan = unique(plan,dims=2)

        println("")
        println("Iteration ", i)
        println("Evaluating function")
        #Evaluate the expensive function
        samples = mapslices(f,plan,dims=1)
        @show minimum(samples)
    end

    
    sol_hist_fitness = Array{Float64,1}()
    sol_hist_iteration = Array{Int64,1}()
    f_c = function (oc)
        push!(sol_hist_fitness, best_fitness(oc))
        push!(sol_hist_iteration, BlackBoxOptim.num_func_evals(oc))
    end        

    res = bboptimize(f; SearchRange=(-5.0,5.0), NumDimensions = 2,
                        PopulationSize=num_start_samples, MaxFuncEvals=10000,
                        CallbackFunction = f_c,CallbackInterval = eps(), 
                        TraceMode = :silent);
    
    equal_iterations_ind = findfirst(x-> x >= length(samples),sol_hist_iteration)
    equal_iterations_fitness = sol_hist_fitness[equal_iterations_ind]

    iterations_ind_for_equal_performance = findfirst(x-> x <= minimum(samples),sol_hist_fitness)
    iterations_for_equal_performance = sol_hist_iteration[iterations_ind_for_equal_performance]

    println("Function evaluations, ", length(samples), ".  Best DE = ", repr(equal_iterations_fitness), ". Best surrogate = ", minimum(samples))
    println("DE iterations for >= performance = ", repr(iterations_for_equal_performance), ".")
    println("Performance increase = ", repr(round(iterations_for_equal_performance/length(samples); digits=2)), ".")
    
    

    return samples, plan, sm_interpolant
end


