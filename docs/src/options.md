# Options
The available options are listed below. The performance using the default options
has been good for several test functions, with and without noise. Make sure to
set `num_interpolants` to a value which is a multiple of the number of processes
used for the least amount of time to create the surrogate. 

```julia
num_start_samples::Int = 5                              # Samples included in the LHC sampling plan
trace::Bool = true                                      # Print progress
sampling_plan_opt_gens::Int = 5000                      # Iterations used to optimize the LHC sampling plan
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
cond_max::Float64 = 1e4                                 # Maximum allowed condition number of RBF matrix A
max_rbf_width::Float64 = 1000.0                         # Maximum RBF width factor
min_rbf_width::Float64 = 1e-4                           # Minimum RBF width factor
max_scale::Float64 = 10.0                               # Maximum linear input dimension scaling factor
min_scale::Float64 = 1e-4                               # Minimum linear input dimension scaling factor
num_interpolants::Int = 20                              # Number of interpolants in ensemble
smooth = :single                                        # Apply smoothing factor, useful for functions with noise
                                                        # false turned off, :single one factor for all points,
                                                        # :variable individual factor for each point, 
                                                        # :single_user user supplied value
max_smooth::Float64 = 1.0                               # Maximum smoothing factor 
smooth_user::Float64 = 0.0                              # User supplied smoothing factor, only applied if smooth = :single_user
iterations::Int64 = 10                                  # Number of infill iterations to be run
num_infill_points::Int64 = 1                            # Number of infill points per iteration
parallel_surrogate::Bool = true                         # Create each surrogate in the ensemble in parallel
infill_funcs::Array{Symbol,1} = [:median,:std]          # Infill criteria, cycled through 
infill_iterations::Int64 = 25_000                       # Iterations to add infill points to the design space
```