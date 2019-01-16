"""
    function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}, options=options())

Optimize the function `f` in the range `search_range` using Radial Basis Function based surrogate model.
"""
function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}, options=options())

    options = _update_options(search_range;options...)

    #Load some optional argument values
    @unpack num_start_samples, sampling_plan_opt_gens,
            iterations, trace = options
    
    #Create sampling plan
    plan = _LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens,trace)
    
    #Evaluate sampling plan
    lhc_samples = f_opt_eval(f,plan,trace)
    criteria = 1

    #Initialize values to be returned
    sm_interpolant = nothing; infill_type = nothing; infill_prediction = nothing
    optres = nothing
    infill_plan = Array{Float64}(undef,size(plan,1),0)
    infill_sample = Array{Float64}(undef,1,0)

    #Run the optimization iterations number of times
    for i = 1:iterations
        if trace
            print("\n \n \n \t Iteration ")
            printstyled(i,bold=true)
            print(" out of ", iterations, "\n")
        end


        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [plan infill_plan]
        sm_interpolant, optres = surrogate_model(samples_all, plan_all, options)
        
        #Points to add to the sampling plan to improve the interpolant
        infill_plan_new, criteria, infill_type_new, infill_prediction_new  = model_infill(plan_all,samples_all,sm_interpolant,criteria,options)
        
        #Evaluate the new infill points
        infill_sample_new = f_opt_eval(f,infill_plan_new,samples_all,trace)

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [infill_sample infill_sample_new]
        infill_type = [infill_type; infill_type_new]
        infill_prediction = [infill_prediction; infill_prediction_new]

    end   
    
    return lhc_samples, plan, sm_interpolant, optres, infill_sample, infill_type, infill_plan, infill_prediction
end


function f_opt_eval(f,plan,samples,trace)

    if trace
        println("Evaluating function ",size(plan,2)," times ...")
    end

    new_samples = mapslices(f,plan,dims=1) 

    if trace
        _, min_loc = findmin(new_samples)
        for i = 1:size(plan,2)
            if i == min_loc[2]
                printstyled(@sprintf("%-15.7g",new_samples[i]); color=:light_green, bold=true)
            else
                printstyled(@sprintf("%-15.7g",new_samples[i]))
            end
        end
        print("\t actual value\n")
        println("---------------------------------------------------------------")

        new_min = minimum(new_samples)
        old_min = minimum(samples)
        new_max = maximum(new_samples)
        old_max = maximum(samples)

        print("Max and min sample value: ")
        printstyled(@sprintf("%.7g",maximum((new_max,old_max))); color=:light_red)
        print("\t")
        printstyled(@sprintf("%.7g",minimum((new_min,old_min))); color=:green, bold=true)

        if isless(new_min,old_min)
            print("\t (Improvement from last iteration ", @sprintf("%.7g",old_min-new_min),")")
        else
            print("\t (Improvement from last iteration N/A)")
        end
        print("\n")
    end

    return new_samples
end

function f_opt_eval(f,plan,trace)

    if trace
        println("Evaluating function ",size(plan,2)," times ...")
    end

    new_samples = mapslices(f,plan,dims=1) 

    if trace
        new_min = minimum(new_samples)
        new_max = maximum(new_samples)


        print("Max and min sample value: ")
        printstyled(@sprintf("%.7g",new_max); color=:light_red)
        print("\t")
        printstyled(@sprintf("%.7g",new_min); color=:light_green, bold=true)
        print("\n")
    end

    return new_samples
end


function compare_de_to_surrogate()
    sol_hist_fitness = Array{Float64,1}()
    sol_hist_iteration = Array{Int64,1}()
    f_c = function (oc)
        push!(sol_hist_fitness, best_fitness(oc))
        push!(sol_hist_iteration, BlackBoxOptim.num_func_evals(oc))
    end        

    res = bboptimize(f; SearchRange=search_range,
                        PopulationSize=num_start_samples, MaxFuncEvals=100000,
                        CallbackFunction = f_c,CallbackInterval = eps(), 
                        TraceMode = :silent);
    
    equal_iterations_ind = findfirst(x-> x >= length(samples),sol_hist_iteration)
    equal_iterations_fitness = sol_hist_fitness[equal_iterations_ind]

    iterations_ind_for_equal_performance = findfirst(x-> x <= minimum(samples),sol_hist_fitness)
    if iterations_ind_for_equal_performance == nothing
        iterations_for_equal_performance = Inf
    else    
        iterations_for_equal_performance = sol_hist_iteration[iterations_ind_for_equal_performance]
    end

    println("Function evaluations, ", length(samples), ".  Best DE = ", repr(equal_iterations_fitness), ". Best surrogate = ", minimum(samples))
    println("DE iterations for >= performance = ", repr(iterations_for_equal_performance), ".")
    println("Performance increase = ", repr(round(iterations_for_equal_performance/length(samples); digits=2)), ".")
end