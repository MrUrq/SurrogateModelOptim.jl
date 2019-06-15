using SurrogateModelOptim
using PlotlyJS
using Statistics
using Distributed

if !(@isdefined loaded)
    loaded = true    
    addprocs(20)
end 
@everywhere using SurrogateModelOptim

# Optimize the test function
Forrester(x) = (6*x[1] - 2)^2 * sin(12*x[1] - 4)
Forrester_disc(x) = x[1]>0.15 ? Forrester(x[1])+15 : Forrester(x[1])
sr=[(0.0,1.0)]

result = smoptimize(Forrester_disc, sr,
                    SurrogateModelOptim.Options(
                    iterations=2,
                    num_interpolants=20, #Preferably even number of added processes
                    num_start_samples=5,
                    rbf_opt_gens=50,
                    infill_iterations=50,
                    num_infill_points=1,
                    trace=true,
                    cond_max=1e10,
                        ));


function plot_fun_1D(fun_original,fun_estimate,sr)    
    N = 101    
    x = range(sr[1][1], length=N,stop=sr[1][2])
    y_original = fun_original.(x)
    y_estimate = [median(fun_estimate([x])) for x in x]
   
    trace1 = scatter(;x=x, y=y_original, mode="lines",name="Original function")
    trace2 = scatter(;x=x, y=y_estimate, mode="lines",name="Estimated function")
    p = plot([trace1,trace2])
end

# Plot the results
display(plot_fun_1D(Forrester_disc,result.sm_interpolant,sr))

return true