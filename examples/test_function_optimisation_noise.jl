using SurrogateModelOptim
using PlotlyJS
using Statistics
using Random

# Optimize the test function Rosenbrock
function rosenbrock_2D(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
search_range=[(-5.0,5.0),(-5.0,5.0)]
fun_max = 90036.0
fun_min = 0.0
noiselevel = (fun_max-fun_min)*0.1 # 10% noise

noisy_rosenbrock_2D = function (x)
    rosenbrock_2D(x)-noiselevel/2+noiselevel*rand()
end

result = smoptimize(noisy_rosenbrock_2D, search_range;
                    options=SurrogateModelOptim.Options(
                    iterations=15,
                    num_interpolants=10, #Preferably even number of added processes
                    num_start_samples=5,
                    smooth=:single,
                        ));

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

display(plot_fun_2D(noisy_rosenbrock_2D,search_range,"Original function with added noise"))
display(plot_fun_2D(rosenbrock_2D,search_range,"Original function without noise"))
display(plot_fun_2D(x->median(result.sm_interpolant(x)),search_range,"Estimated function"))