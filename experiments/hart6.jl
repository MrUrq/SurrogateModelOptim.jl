using PGFPlotsX
using Colors
using ColorSchemes
using Random
using SurrogateModelOptim
using ColorBrewer
using Distances
using ScatteredInterpolation
using Statistics

res = smoptimize(hart_6D, repeat([(0.0,1.0)],6),
               SurrogateModelOptim.options(num_interpolants=20,
                                           num_start_samples=10,
                                           sampling_plan_opt_gens=100_000,
                                           rbf_opt_gens=50_000,
                                           variable_kernel_width = true,
                                           variable_dim_scaling = true,
                                           smooth = :single,
                                           num_infill_points = 1,
                                           infill_funcs = [:median,:confint,:median,:std,:median,:dist],
                                           iterations =150,
                                           trace = true))
return res;



