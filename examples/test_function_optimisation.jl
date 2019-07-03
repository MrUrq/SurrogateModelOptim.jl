using SurrogateModelOptim
using PlotlyJS
using Statistics

# Optimize the test function Rosenbrock
function rosenbrock_2D(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
search_range=[(-5.0,5.0),(-5.0,5.0)]

result = smoptimize(rosenbrock_2D, search_range;
                    options=SurrogateModelOptim.Options(
                    iterations=15,
                    num_interpolants=10, #Preferably even number of added processes
                    num_start_samples=5,
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

display(plot_fun_2D(rosenbrock_2D,search_range,"Original function"))
display(plot_fun_2D(x->median(result.sm_interpolant(x)),search_range,"Estimated function"))