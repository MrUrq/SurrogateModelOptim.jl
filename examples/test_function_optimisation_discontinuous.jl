using SurrogateModelOptim
using PlotlyJS
using Statistics

# Discontinuous Forrester function
Forrester(x) = (6*x[1] - 2)^2 * sin(12*x[1] - 4)
Forrester_disc(x) = x[1]>0.15 ? Forrester(x[1])+15 : Forrester(x[1])
search_range=[(0.0,1.0)]

# Optimize the test function
result = smoptimize(Forrester_disc, search_range;
                    options=SurrogateModelOptim.Options(
                    iterations=20,
                    num_interpolants=10, #Preferably even number of added processes
                    num_start_samples=5,
                        ));

function plot_fun_1D(fun_original,fun_estimate,search_range)    
    N = 101    
    x = range(search_range[1][1], length=N,stop=search_range[1][2])
    y_original = fun_original.(x)
    y_estimate = [median(fun_estimate([x])) for x in x]
   
    trace1 = scatter(;x=x, y=y_original, mode="lines",name="Original function")
    trace2 = scatter(;x=x, y=y_estimate, mode="lines",name="Estimated function")
    p = plot([trace1,trace2])
end

# Plot the results
display(plot_fun_1D(Forrester_disc,result.sm_interpolant,search_range))