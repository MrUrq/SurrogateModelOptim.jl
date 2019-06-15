using SurrogateModelOptim
using PlotlyJS
using Statistics
dir_path = @__DIR__ 
include(joinpath(dir_path,"test_functions.jl"))

# Optimize the test function
func = test_funs[:rosenbrock_2D]
result = smoptimize(func.fun, func.sr,
                    SurrogateModelOptim.Options(
                    iterations=7,
                    num_interpolants=2, #Preferably even number of added processes
                    num_start_samples=5,
                    rbf_opt_gens=50_000,
                    infill_iterations=50_000,
                    num_infill_points=3,
                    trace=true,
                        ));

# This runs num_start_samples + (iterations*num_infill_points) function
# evaluations in total.

function plot_fun_2D(fun,sr,title)    
    N = 51    
    x = range(sr[1][1], stop = sr[1][2], length = N)
    y = range(sr[2][1], stop = sr[2][2], length = N)

    grid(x,y) = [x,y]
    z = @. fun(grid(x,y'))

    trace = surface(x=x,y=y,z=z, colorscale="Viridis")     
    layout = Layout(title=title)
    p = plot(trace,layout) 
end

# Plot the results if the optimised function is 2-dimensional
if length(func.sr) == 2
    display(plot_fun_2D(func.fun,func.sr,"Original function"))
    display(plot_fun_2D(x->median(result.sm_interpolant(x)),func.sr,"Estimated function"))
end

return true