"""
    smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}; options=Options())

Optimize the function `f` in the range `search_range` using a Radial Basis Function based surrogate model.
"""
function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}; options::Options=Options())

    #Load some option values
    @unpack num_start_samples, sampling_plan_opt_gens,
            iterations, trace, create_final_surrogate = options
    
    #Create sampling plan
    lhc_plan = scaled_LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens;trace=trace)
    
    #Evaluate sampling plan
    lhc_samples = f_opt_eval(f,lhc_plan;trace=trace)

    #Initialize variables to be returned
    sm_interpolant = nothing
    infill_type = Array{Symbol,1}(undef,0)
    infill_prediction = Array{Float64,1}(undef,0)
    optres = nothing
    infill_plan = Array{Float64,2}(undef,size(lhc_plan,1),0)
    infill_sample = Array{Float64,2}(undef,1,0)

    #Run the entire optimization iterations number of times
    for i = 1:iterations
        
        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = surrogate_model(plan_all, samples_all; options=options)
        
        #Points to add to the sampling plan to improve the interpolant
        infill_plan_new, infill_type_new, infill_prediction_new, options  = model_infill(search_range,plan_all,
                                                                                samples_all,sm_interpolant;options=options)
        
        #Evaluate the new infill points
        print_iteration(trace,i,iterations)
        infill_sample_new = f_opt_eval(f,infill_plan_new,samples_all;trace=trace)
        

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [infill_sample infill_sample_new]
        infill_type = [infill_type; infill_type_new]
        infill_prediction = [infill_prediction; infill_prediction_new]

    end

    if create_final_surrogate
        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = surrogate_model(plan_all, samples_all; options=options)
    end
    
    return SurrogateResult( lhc_samples, lhc_plan, sm_interpolant,
                            optres, infill_sample, infill_type,
                            infill_plan, infill_prediction,options)
end

"""
    f_opt_eval(f,plan,samples;trace::Symbol=:silent)

Calculate the objective function value(s) and provide intermediate results printing
showing improvements over previous best iteration.
"""
function f_opt_eval(f,plan,samples;trace::Symbol=:silent)

    (trace == :verbose) && println("Evaluating function ",size(plan,2)," times")

    new_samples = mapslices(f,plan,dims=1) 

    print_f_opt(trace,new_samples,samples,plan)

    return new_samples
end

function print_f_opt(trace,new_samples,samples,plan)

    new_min = minimum(new_samples)
    old_min = minimum(samples)
    new_max = maximum(new_samples)
    old_max = maximum(samples)

    if trace == :compact
        min_sample = minimum((old_min,new_min))
        if isless(new_min,old_min)
            print(@sprintf("%20.7g %20.7g\n", min_sample, old_min-new_min))
        else
            print(@sprintf("%20.7g %20.7s\n", min_sample, "N/A"))
        end
        
    elseif trace == :verbose
        print("    Minimum sample value ")
        printstyled(@sprintf("%.7g",new_min); color=:light_green, bold=true)
        print("\t(maximum = ")
        printstyled(@sprintf("%.7g",new_max);)
        print(")")

        if isless(new_min,old_min)
            print("\t (Improvement from last iteration ", @sprintf("%.7g",old_min-new_min),")\n")
        else
            print("\t (Improvement from last iteration N/A)\n")
        end
        print("\n")
    end
end

"""
    f_opt_eval(f,plan;trace::Symbol=:silent)

Calculate the objective function value(s) and plot value if trace.
"""
function f_opt_eval(f,plan;trace::Symbol=:silent)

    (trace == :verbose) && println("Evaluating function ",size(plan,2)," times")
    (trace == :compact) && print(@sprintf("%s %21s %20s\n", "   Iteration", "Function value", "Improvement"))

    new_samples = mapslices(f,plan,dims=1) 

    if trace == :verbose
        new_min = minimum(new_samples)
        new_max = maximum(new_samples)

        print("    Minimum sample value ")
        printstyled(@sprintf("%.7g",new_min); color=:light_green, bold=true)
        print("\t(maximum = ")
        printstyled(@sprintf("%.7g",new_max);)
        println(")\n")
    end

    return new_samples
end

"""
    print_iteration(i,iterations)

REPL printing of the current iteration if trace.
"""
function print_iteration(trace,i,iterations)
    if trace == :compact
        print(@sprintf("%4.4g%s%g",i, " out of ", iterations))
    end

    if trace == :verbose
        print("\n\nFunction eval. ")
        printstyled("iteration ";bold=true,color=:light_red)
        printstyled(i,bold=true)
        print(" out of ", iterations)
        if trace == :verbose
            println("\n")
        end
    end
end