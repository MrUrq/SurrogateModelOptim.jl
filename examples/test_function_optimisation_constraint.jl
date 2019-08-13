using LatinHypercubeSampling
using PlotlyJS
using Statistics
using SurrogateModelOptim
using Parameters

options=SurrogateModelOptim.Options(
    iterations=10,
    num_interpolants=10, #Preferably even number of added processes
    num_start_samples=5,
    create_final_surrogate=true, #Use the results from last iteration to
                                 #re-create the surrogate before using it for plotting
        )
    
function rosenbrock_2D(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

rosenbrock_2D_pen = function (x)            
    if abs(x[1]) < 1
        return rosenbrock_2D(x).+50_000 #Add large value as penalty
    else
        return rosenbrock_2D(x)
    end
end

search_range=[(-5.0,5.0),(-5.0,5.0)]

function constrained_smoptimize(f,search_range,options)

    #Load some option values
    @unpack num_start_samples, sampling_plan_opt_gens,
            iterations, trace = options
    
    #Create sampling plan
    lhc_plan = SurrogateModelOptim.scaled_LHC_sampling_plan(search_range,num_start_samples,10000)
    
    #Evaluate sampling plan
    lhc_samples = SurrogateModelOptim.f_opt_eval(f,lhc_plan;trace=trace)

    #Initialize variables to be returned
    sm_interpolant_pen = nothing
    infill_type = Array{Symbol,1}(undef,0)
    infill_prediction = Array{Float64,1}(undef,0)
    optres = nothing
    infill_plan = Array{Float64,2}(undef,size(lhc_plan,1),0)
    infill_sample = Array{Float64,2}(undef,1,0)

    #Run the entire optimization iterations number of times
    for i = 1:iterations
        SurrogateModelOptim.print_iteration(trace,i,iterations)
        
        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = SurrogateModelOptim.surrogate_model(plan_all, samples_all; options=options)

        #Add the penalty if abs(x) < 1
        sm_interpolant_pen = function (x)            
            if abs(x[1]) < 1
                return sm_interpolant(x).+50_000 #Add large value as penalty
            else
                return sm_interpolant(x)
            end
        end
        
        #Points to add to the sampling plan to improve the interpolant
        infill_plan_new, infill_type_new, infill_prediction_new, options  = 
            model_infill(search_range,plan_all,samples_all,sm_interpolant_pen;options=options)
            
        #Evaluate the new infill points
        infill_sample_new = SurrogateModelOptim.f_opt_eval(f,infill_plan_new,samples_all;trace=trace)

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [infill_sample infill_sample_new]
        infill_type = [infill_type; infill_type_new]
        infill_prediction = [infill_prediction; infill_prediction_new]
    end   
    
    return SurrogateModelOptim.SurrogateResult( lhc_samples, lhc_plan, sm_interpolant_pen,
                            optres, infill_sample, infill_type,
                            infill_plan, infill_prediction,options)
end
result = constrained_smoptimize(rosenbrock_2D,search_range,options)
show(result)

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

# Plot the results 
display(plot_fun_2D(rosenbrock_2D_pen,search_range,"Original function with added penalty"))
display(plot_fun_2D(x->median(result.sm_interpolant(x)),search_range,"Estimated function with added penalty"))