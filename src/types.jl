@with_kw struct Options
    num_start_samples::Int = 10
    show_trace::Bool = true
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
    max_scale::Float64 = 10.0
    min_scale::Float64 = 1e-4
    num_interpolants::Int = 1
    smooth = false
    max_smooth = 1.0
end


"""
    RBFHypers(width::S,kernelFunc,scaling::U,fitness::Float64)

Datastructure to store results from the optimisation of an RBF interpolation kernel
"""
struct RBFHypers{T,U}
    width::T
    kernelFunc
    scaling::U
    fitness::Float64
end
# struct RBFHypers{U,V}
#     kernelFunc
#     scaling::U
#     smooth::V
# end

struct SurrogateEstimate{T}
    sm_estimate::T
end