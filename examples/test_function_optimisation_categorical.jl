using LatinHypercubeSampling
using PlotlyJS
using Statistics
using SurrogateModelOptim
using Parameters

options=SurrogateModelOptim.Options(
    iterations=15,
    num_interpolants=10, #Preferably even number of added processes
    num_start_samples=5,
        )
    
function rosenbrock_2D(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

search_range=[(-5.0,5.0),(-5.0,5.0)]

possible_locs = 20 #Defined the possible design locations from -5,5
dims = [LatinHypercubeSampling.Categorical(possible_locs) for i in 1:length(search_range)]
vals = [Tuple((i for i in range(-5.0,stop=5.0,length=possible_locs))) for i in 1:length(search_range)] 

function categorical_smoptimize(f,search_range,dims,vals,options)

    #Print the known minimum (brute force approach)
    possible_designs = collect(Base.Iterators.product(vals...))    
    println("Number of possible designs = ", length(possible_designs))
    min_loc = argmin(f.(possible_designs))
    max_loc = argmax(f.(possible_designs))
    println("Minimum and maximum value in categorical design space:")
    println("Min ",f(possible_designs[min_loc]), ", ", possible_designs[min_loc])
    println("Max ",f(possible_designs[max_loc]), ", ", possible_designs[max_loc])

    #Load some option values
    @unpack num_start_samples, sampling_plan_opt_gens,
            iterations, trace = options
    
    #Create categorical sampling plan with the values in vals
    initialSample = randomLHC(num_start_samples,dims)
    Xind = permutedims(LHCoptim!(initialSample,sampling_plan_opt_gens;dims=dims)[1])
    lhc_plan = Array{Float64,2}(undef,size(Xind)) 
    for i = 1:size(Xind,1)
        for j = 1:size(Xind,2)
            lhc_plan[i,j] = vals[i][Xind[i,j]]
        end
    end
    
    #Evaluate sampling plan
    lhc_samples = SurrogateModelOptim.f_opt_eval(f,lhc_plan;trace=trace)

    #Initialize variables to be returned
    sm_interpolant = nothing
    infill_type = Array{Symbol,1}(undef,0)
    infill_prediction = Array{Float64,1}(undef,0)
    optres = nothing
    infill_plan = Array{Float64,2}(undef,size(lhc_plan,1),0)
    infill_sample = Array{Float64,2}(undef,1,0)

    #Run the entire optimization iterations number of times
    for i = 1:iterations
        trace && SurrogateModelOptim.print_iteration(i,iterations)
        
        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = SurrogateModelOptim.surrogate_model(plan_all, samples_all; options=options)

        #Choose the closest categorical value when evaluating the function. Note
        # this alters the input of x to match closest value.
        sm_interpolant_cat! = function (x,vals)            
            for (i,x_val) in enumerate(x)
                c_ind = SurrogateModelOptim.closest_index(x_val, vals[i])
                x[i] = vals[i][c_ind]
            end
            return sm_interpolant(x)
        end
        
        #Points to add to the sampling plan to improve the interpolant
        infill_plan_new, infill_type_new, infill_prediction_new, options  = 
            model_infill(search_range,plan_all,samples_all,x->sm_interpolant_cat!(x,vals);options=options)
            
        #Evaluate the new infill points
        infill_sample_new = SurrogateModelOptim.f_opt_eval(f,infill_plan_new,samples_all;trace=trace)

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [infill_sample infill_sample_new]
        infill_type = [infill_type; infill_type_new]
        infill_prediction = [infill_prediction; infill_prediction_new]
    end   
    
    return SurrogateModelOptim.SurrogateResult( lhc_samples, lhc_plan, sm_interpolant,
                            optres, infill_sample, infill_type,
                            infill_plan, infill_prediction,options)
end
result = categorical_smoptimize(rosenbrock_2D,search_range,dims,vals,options)

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
display(plot_fun_2D(rosenbrock_2D,search_range,"Original function"))
display(plot_fun_2D(x->median(result.sm_interpolant(x)),search_range,"Estimated function"))