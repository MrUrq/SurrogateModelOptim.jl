# function for creating optimised surrogate
function smoptimize(f::Function; SearchRange, NumDimensions=false, kwargs...)

    params = _setup_parameters(SearchRange,NumDimensions; kwargs...)

    #Load the optional argument values
    @Parameters.unpack NumStartSamples, TraceMode, SamplingPlanOptGens = params
    
    #Create sampling plan
    plan = _LHC_sampling_plan(SearchRange,NumDimensions,NumStartSamples,SamplingPlanOptGens)

    #

end




     

#     #evaluate f in the LHC points

    
#     #surrogate model creation
#     smsetup(samples, points; kwargs...)

#     #infill points
#     surrogate_infill(sm,samples,points)

# end

# function for creating an optimised surrogate
function smsetup(samples,points;)
    # 
end




function testf(;kwargs...)
    
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
# Default parameters for all convenience methods that are exported to the end user.
# See `OptRunController` for the description.
# """
# const DefaultParameters = ParamsDict(
#     :NumDimensions  => :NotSpecified, # Dimension of problem to be optimized
#     :SearchRange    => (-1.0, 1.0), # Default search range to use per dimension unless specified
#     :SearchSpace    => false, # Search space can be directly specified and will then take precedence over NumDimensions and SearchRange.
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

#     :TraceMode      => :compact,  # Print tracing information during the optimization
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