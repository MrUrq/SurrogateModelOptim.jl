using PGFPlotsX
using Colors
using ColorSchemes
using Random
using SurrogateModelOptim
using ColorBrewer
using Distances
using ScatteredInterpolation
using Statistics

res = smoptimize(hart6, repeat([(0.0,1.0)],6),
               SurrogateModelOptim.options(num_interpolants=20,
                                           num_start_samples=5,
                                           sampling_plan_opt_gens=10000,
                                           rbf_opt_gens=10000,
                                           variable_kernel_width = true,
                                           variable_dim_scaling = true,
                                           smooth = :single,
                                           smooth_user = 0.0))
return nothing       



