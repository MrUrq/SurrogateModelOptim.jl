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
                                           num_start_samples=15,
                                           sampling_plan_opt_gens=25_000,
                                           rbf_opt_gens=50_000,
                                           variable_kernel_width = true,
                                           variable_dim_scaling = true,
                                           smooth = :single,
                                           num_infill_points = 1,
                                           infill_funcs = [:min,:min_2std,:std,:dist],
                                           iterations = 10,
                                           trace = true))
return res;



