using SurrogateModelOptim
using PlotlyJS
using Statistics
using Random
dir_path = @__DIR__ 
include(joinpath(dir_path,"test_functions.jl"))

func = test_funs[:rosenbrock_2D]
noiselevel = (func.max_val-func.min_val)*0.1

# Add deterministic noise which gives the same value each time the 
# same point is sampled.
opt_fun = function (x)
    noise = 0
    try
        noise = noiselevel*rand(MersenneTwister(abs(sum(reinterpret(Int64,x)))))
    catch 
        noise = noiselevel*rand(MersenneTwister(1))
    end
    func.fun(x)-noiselevel/2+noise
end

# Optimize the test function
result = smoptimize(opt_fun, func.sr;
                    options=SurrogateModelOptim.Options(
                    iterations=5,
                    num_interpolants=2, #Preferably even number of added processes
                    num_start_samples=5,
                    rbf_opt_gens=50,
                    infill_iterations=50,
                    num_infill_points=2,
                    trace=true,
                        ));

# This runs num_start_samples + (iterations*num_infill_points) function
# evaluations in total.


function plot_fun_2D(fun,sr,title)    
    N = 100    
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
    display(plot_fun_2D(opt_fun,func.sr,"Original function"))
    display(plot_fun_2D(x->median(result.sm_interpolant(x)),func.sr,"Estimated function"))
end

return true