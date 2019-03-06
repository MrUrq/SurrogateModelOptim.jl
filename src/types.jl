"""
    Options

Options configurable by the user with recommended default values.
"""
@with_kw struct Options
    num_start_samples::Int = 10
    trace::Bool = true
    sampling_plan_opt_gens::Int = 1000
    rippa::Bool = true
    kerns = [ScatteredInterpolation.Gaussian, ScatteredInterpolation.InverseQuadratic,
             ScatteredInterpolation.InverseMultiquadratic]
    rbf_opt_gens::Int = 25_000
    rbf_opt_method::Symbol = :de_rand_1_bin_radiuslimited
    rbf_dist_metric = Distances.Euclidean()
    variable_kernel_width::Bool = true
    variable_dim_scaling::Bool = true
    cond_max::Float64 = 1e4
    max_rbf_width::Float64 = 1000.0
    min_rbf_width::Float64 = 1e-4
    max_scale::Float64 = 10.0
    min_scale::Float64 = 1e-4
    num_interpolants::Int = 20
    smooth = :single
    max_smooth::Float64 = 1.0
    smooth_user::Float64 = 0.0
    iterations::Int64 = 10
    num_infill_points::Int64 = 1
    parallel_surrogate::Bool = true
    infill_funcs::Array{Symbol,1} = [:std,:median]
    infill_iterations::Int64 = 10_000
    categorical::Bool = false

    @assert ((smooth == false) || (smooth == :variable) || 
    (smooth == :single) || (smooth == :single_user)) "Not supported option for smooth"
end


"""
    RBFHypers(width::S,kernelFunc,scaling::U,fitness::Float64)

Data structure to store results from the optimisation of an RBF interpolation kernel
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




