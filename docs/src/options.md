# Options
The available options are listed below. The performance using the default options has been
good for several test functions, with and without noise. Make sure to set
`num_interpolants` to a value which is a multiple of the number of processes used for the
least amount of time to create the surrogate. Based on a small benchmark test, the number
of interpolants has marginal improvement to the optimisation above 10 interpolants. 

```julia
num_start_samples::Int = 5                              # Samples included in the LHC sampling plan
trace::Bool = true                                      # Print progress
sampling_plan_opt_gens::Int = 50_000                    # Iterations used to optimize the LHC sampling plan
rippa::Bool = true                                      # Rippas algorithm to reduce computational effort optimising the surrogate
kerns = [ScatteredInterpolation.Gaussian,               # RBF kernels to choose from
        ScatteredInterpolation.InverseQuadratic,        
        ScatteredInterpolation.InverseMultiquadratic]
rbf_opt_gens::Int = 50_000                              # Generations that the RBF hyperparameters are optimized
rbf_opt_pop::Int = 50                                   # Population size of RBF hyperparameter optimization
rbf_opt_method::Symbol = :de_rand_1_bin_radiuslimited   # BlackBoxOptim optimization method for RBF hyperparameters
rbf_dist_metric = Distances.Euclidean()                 # Distance metric used to create the RBF
variable_kernel_width::Bool = true                      # Allow individual kernels and widths for each sample
variable_dim_scaling::Bool = true                       # Linearly scale the input dimensions for best fit
cond_max::Float64 = 1e10                                # Maximum allowed condition number of RBF matrix A
max_rbf_width::Float64 = 1.0                            # Maximum RBF width factor
min_rbf_width::Float64 = 1e-10                          # Minimum RBF width factor
max_scale::Float64 = 1.0                                # Maximum linear input dimension scaling factor
min_scale::Float64 = 1e-10                              # Minimum linear input dimension scaling factor
num_interpolants::Int = 10                              # Number of interpolants in ensemble
smooth = false                                          # Apply smoothing factor, useful for functions with noise
                                                        # false turned off, :single one factor for all points,
                                                        # :variable individual factor for each point, 
                                                        # :single_user user supplied value
max_smooth::Float64 = 0.005                             # Maximum smoothing factor 
smooth_user::Float64 = 0.0                              # User supplied smoothing factor, only applied if smooth = :single_user
iterations::Int64 = 5                                   # Number of infill iterations to be run
num_infill_points::Int64 = 1                            # Number of infill points per iteration
parallel_surrogate::Bool = true                         # Create each surrogate in the ensemble in parallel
infill_funcs::Array{Symbol,1} = [:median,:std]          # Infill criteria, cycled through 
infill_iterations::Int64 = 50_000                       # Iterations to add infill points to the design space
```


## Smoothing
`smooth`, `max_smooth` and `smooth_user` are all related to the ability of smoothing the
RBF interpolant to handle cases with noise. The best performance is had when the smallest
amount of smoothing is applied while being large enough to hinder high frequency
oscillations in the resulting interpolant. If the input is noisy, set `smooth` to
`:single` and the amount of smoothing applied is automatically optimised. The amount of
smoothing is bound between 0 and `max_smooth` when optimised automatically. It can also be
set to a fixed value with the option `smooth = :single_user` where the value is supplied
with `smooth_user` option.

## Kernel options
`variable_kernel_width` is by default `true` meaning that the RBF kernel and the width of
each kernel is optimised individually for each sampled point. This does increase the
degree of freedom which can lead to overfitting when the number of sample points is small.
The performance when using `variable_kernel_width = true` is often smaller initially but
can converged the problem to lower values. Recommendation is to leave the option `true` at
all times.

Radial Basis Function interpolation is as the name suggests, radial in space. This can
create problems when vastly differing input to output dimensions are used such as in the
2D Rosenbrock benchmark problem. In such cases the optimal radial width becomes a
compromise between the two input dimensions. To solve this, the input dimensions can be
linearly scaled automatically using `variable_dim_scaling = true`, which is the default. 
