"""
    Options

Options configurable by the user with recommended default values.
"""
@with_kw struct Options
    num_start_samples::Int = 5
    trace::Symbol = :compact
    sampling_plan_opt_gens::Int = 1_000
    rippa::Bool = true
    kerns = [ScatteredInterpolation.Gaussian]
    rbf_opt_gens::Int = 10_000
    rbf_opt_pop::Int = 50
    rbf_opt_method::Symbol = :de_rand_1_bin_radiuslimited
    rbf_dist_metric = Distances.Euclidean()
    variable_kernel_width::Bool = true
    variable_dim_scaling::Bool = true
    cond_max::Float64 = 1e8
    max_rbf_width::Float64 = 1.0
    min_rbf_width::Float64 = 1e-10
    max_scale::Float64 = 1.0
    min_scale::Float64 = 1e-10
    num_interpolants::Int = 10
    smooth = false
    max_smooth::Float64 = 0.005
    smooth_user::Float64 = 0.0
    iterations::Int64 = 5
    num_infill_points::Int64 = 1
    parallel_surrogate::Bool = true
    infill_funcs::Array{Symbol,1} = [:median,:std]
    infill_iterations::Int64 = 10_000
    create_final_surrogate::Bool = false

    @assert ((smooth == false) || (smooth == :variable) || 
    (smooth == :single) || (smooth == :single_user)) "Not supported smoothing option"

    @assert ((trace == :silent) || (trace == :compact) || 
    (trace == :verbose)) "Not supported trace option"
end


"""
    RBFHypers(width::S,kernelFunc,scaling::U,fitness::Float64)

Data structure to store results from the optimization of an RBF interpolation kernel
"""
struct RBFHypers{T,U}
    kernelFunc
    scaling::T
    smooth::U
end


"""
    SurrogateResult

Data structure to store the results of an optimization task. 
"""
struct SurrogateResult
    lhc_samples::Array{Float64,2}
    lhc_plan::Array{Float64,2}
    sm_interpolant
    sm_interpolant_settings::Array{SurrogateModelOptim.RBFHypers{T,U},1} where U where T
    infill_samples::Array{Float64,2}
    infill_type::Array{Symbol,1}
    infill_plan::Array{Float64,2}
    infill_prediction::Array{Float64,1}
    options::Options
end

function Base.show(io::IO, res::SurrogateResult)
    println("Surrogate Model Optim Result")
    println(io, @sprintf("  Best fitness: \t%5.7g,  (worst %5.7g)",best_fitness(res),maximum([res.lhc_samples res.infill_samples])))
    println(io, "  Best candidate: \t",best_candidate(res))    
    println(io, "  Function calls: \t",f_calls(res))
    println(io, "  Iterations: \t\t",res.options.iterations)
    println(io, "  LHC sampling points: \t",size(res.lhc_plan,2))
    println(io, "  Infill criteria: \t",unique(res.infill_type))
    if res.sm_interpolant_settings[1].scaling != false
        sc = mean([opt.scaling./opt.scaling[1] for opt in res.sm_interpolant_settings])
        println(io, "  Mean axis scaling: \t",sc," (relative to first dimension)")
    else
        println(io, "  Mean axis scaling: \t","scaling is fixed")
    end
    println(io, "  Smooth: \t\t",res.options.smooth, " (set to :single if noise is expected)")
    print(io, "  Returned surrogate contains all samples: ",res.options.create_final_surrogate)
end