@with_kw struct Options
    num_start_samples::Int = 10
    trace::Bool = true
    sampling_plan_opt_gens::Int = 1000
    rippa::Bool = true
    kerns = [ScatteredInterpolation.Gaussian, ScatteredInterpolation.InverseQuadratic,
             ScatteredInterpolation.InverseMultiquadratic]
    rbf_opt_gens::Int = 1000
    rbf_opt_method::Symbol = :de_rand_1_bin_radiuslimited
    rbf_dist_metric = Distances.Euclidean()
    variable_kernel_width::Bool = false
    variable_dim_scaling::Bool = true
    cond_max::Float64 = 1e4
    max_rbf_width::Float64 = 1000.0
    min_rbf_width::Float64 = 1e-4
    max_scale::Float64 = 10.0
    min_scale::Float64 = 1e-4
    num_interpolants::Int = 1
    smooth = false
    max_smooth::Float64 = 0.01
    smooth_user::Float64 = 0.0
    iterations::Int64 = 10
    num_infill_points::Int64 = 1
    infill_funcs::Array{Symbol,1} = [:min,:min_2std,:std]
end


"""
    RBFHypers(width::S,kernelFunc,scaling::U,fitness::Float64)

Datastructure to store results from the optimisation of an RBF interpolation kernel
"""
struct RBFHypers{T,U}
    kernelFunc
    scaling::T
    smooth::U
end


struct SurrogateEstimate{T}
    sm_estimate::T
end