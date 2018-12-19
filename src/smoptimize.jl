# function for creating optimised surrogate
function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}, options=options())

    options = _update_options(search_range;options...)

    #Load the optional argument values
    @unpack num_start_samples, show_trace, sampling_plan_opt_gens = options
    
    #Create sampling plan
    plan = _LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens)

    #Evaluate the expensive function
    samples = mapslices(f,plan,dims=1)

    #Create the optimized Radial Basis Function interpolant 
    rbf_hypers = surrogate_model(samples, plan, options)
    
    

    # #     #infill points
    # #     surrogate_infill(sm,samples,points)
    # # Loop over everything
end




# function for creating an optimised surrogate
function smsetup(samples,points;)
    # 
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



