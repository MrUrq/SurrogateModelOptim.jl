using PGFPlotsX
using Colors
using ColorSchemes
using Random
using SurrogateModelOptim
using ColorBrewer
using Distances
using ScatteredInterpolation
using Statistics


f = rosenbrock2d
noiselevel = 10e4*0.06
#rotatedHyperElipsoid
#styblinskiTang

res = smoptimize(f, [(-5.0,5.0),(-5.0,5.0)],
        SurrogateModelOptim.options(num_interpolants=20,
                                    num_start_samples=20,
                                    sampling_plan_opt_gens=100000,
                                    rbf_opt_gens=10000,
                                    variable_kernel_width = true,
                                    variable_dim_scaling = true,
                                    smooth = :single,
                                    smooth_user = 0.0));display(funcplotX(x->median(res[3](x))[1],-5,5,-5,5))



res = smoptimize(x->(f(x)-noiselevel/2+noiselevel*rand(MersenneTwister(abs(sum(reinterpret(Int64,x)))))), [(-5.0,5.0),(-5.0,5.0)],
        SurrogateModelOptim.options(num_interpolants=20,
                                    num_start_samples=20,
                                    sampling_plan_opt_gens=100000,
                                    rbf_opt_gens=10000,
                                    variable_kernel_width = true,
                                    variable_dim_scaling = true,
                                    smooth = :single,
                                    smooth_user = 0.0));display(funcplotX(x->median(res[3](x))[1],-5,5,-5,5))
                           
return nothing

