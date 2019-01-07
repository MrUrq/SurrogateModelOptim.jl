

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
    for i = 2:10

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

    # tmp = Array{Float64,1}(10000)
    # f_s = x->push!(tmp,f(best_candidate(x)))
    # res = bboptimize(f;
    #         Method=:de_rand_1_bin,SearchRange=search_range, PopulationSize=num_start_samples, MaxFuncEvals=10000,
    #         CallbackFunction = f_s,CallbackInterval = eps()*0.5),tmp;
    # @show best_de_solution = f(res.archive_output.best_candidate)
    res = bboptimize(f;
            Method=:de_rand_1_bin,SearchRange=search_range, PopulationSize=num_start_samples, MaxFuncEvals=length(samples),
            TraceMode=:silent);
    @show best_de_solution = f(res.archive_output.best_candidate)


    ftol = 20000
    res = bboptimize(f;
            Method=:de_rand_1_bin,SearchRange=search_range, PopulationSize=num_start_samples, TargetFitness=minimum(samples)-ftol,
            FitnessTolerance = ftol, TraceMode=:silent);
    @show equal_de_solution = res.f_calls


            println("Best de = ", best_de_solution, ". Best surrogate = ", minimum(samples), ". Function evals. ", length(samples), ". Iterations to get equally good de solution = ", equal_de_solution, ". Performance inrease = ", equal_de_solution/length(samples))

    return samples, plan, sm_interpolant, best_de_solution, minimum(samples), res
end




# # trace current optimization state,
# # Called by OptRunController trace_progress()
# function trace_state(io::IO, alg::BorgMOEA, mode::Symbol)
#     println(io, "pop.size=", popsize(alg.population),
#                 " arch.size=", length(archive(alg)),
#                 " n.restarts=", alg.n_restarts)
#     if mode == :verbose
#         # output recombination operator rates
#         println(io, "P(recombine):")
#         for i in eachindex(alg.recombinate)
#             println(io, "  #$i(", alg.recombinate[i], ")=",
#                     @sprintf("%.4f",alg.recombinate_distr.p[i]))
#         end
#     end
# end



