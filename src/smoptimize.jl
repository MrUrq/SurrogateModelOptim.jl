# function for creating optimised surrogate
function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}, options=options())

    options = _update_options(search_range;options...)

    #Load the optional argument values
    @unpack num_start_samples, show_trace, sampling_plan_opt_gens = options
    
    #Create sampling plan
    plan = _LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens)

    #Evaluate the expensive function
    samples = mapslices(f,plan,dims=1)

    #Optimize the Radial Basis Function interpolation hyperparameters 
    rbf_hypers = rbf_hypers_opt(samples, plan, options)

    # #     #infill points
    # #     surrogate_infill(sm,samples,points)
    # # Loop over everything
end




# function for creating an optimised surrogate
function smsetup(samples,points;)
    # 
end




function testf(;kwargs...)
    return kwargs
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






# """    


# end
# Default parameters for all convenience methods that are exported to the end user.
# See `OptRunController` for the description.
# """
# const DefaultParameters = ParamsDict(
#     :num_dimensions  => :NotSpecified, # Dimension of problem to be optimized
#     :search_range    => (-1.0, 1.0), # Default search range to use per dimension unless specified
#     :SearchSpace    => false, # Search space can be directly specified and will then take precedence over num_dimensions and search_range.
#     :FitnessScheme  => MinimizingFitnessScheme, # fitness scheme to be used
#     :TargetFitness => nothing, # optimal (target) fitness, if known

#     :Method => :adaptive_de_rand_1_bin_radiuslimited,

#     :MaxTime        => 0.0,
#     :MaxFuncEvals   => 0,
#     :MaxSteps       => 10000,
#     :MaxStepsWithoutProgress => 10000,
#     :MinDeltaFitnessTolerance => 1e-50,
#     :FitnessTolerance => 1e-8,

#     :MaxNumStepsWithoutFuncEvals => 100,

#     :NumRepetitions => 1,     # Number of repetitions to run for each optimizer for each problem

#     :show_trace      => :compact,  # Print tracing information during the optimization
#     :TraceInterval  => 0.50,  # Minimum number of seconds between consecutive trace messages printed to STDOUT
#     :SaveTrace      => false,
#     :SaveFitnessTraceToCsv => false, # Save a csv file with information about the major fitness improvement events (only the first event in each fitness magnitude class is saved)
#     :SaveParameters => false, # Save parameters to a json file for later scrutiny

#     :CallbackFunction => x -> x, # Function to callback to, here just the identity function.
#     :CallbackInterval  => 0.0,  # Minimum number of seconds between consecutive callbacks. If 0.0 we never callback.

#     :RandomizeRngSeed => true, # Randomize the RngSeed value before using any random numbers.
#     :RngSeed        => 1234,   # The specific random seed to set before any random numbers are generated. The seed is randomly selected if RandomizeRngSeed is true, and this parameter is updated with its actual value.

#     :PopulationSize => 50
# )