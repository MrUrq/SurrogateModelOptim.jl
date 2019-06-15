# using Distributed
using LatinHypercubeSampling
using PlotlyJS
using Statistics

# if !(@isdefined loaded)
#     loaded = true    
#     addprocs(20)
# end 
# @everywhere using SurrogateModelOptim
    


dir_path = @__DIR__ 
include(joinpath(dir_path,"test_functions.jl"))

options = SurrogateModelOptim.Options(
    iterations=5, num_interpolants=2, #Preferably even number of added processes
    num_start_samples=4, rbf_opt_gens=50, infill_iterations=50,
    num_infill_points=1, trace=true, categorical=true,
    variable_kernel_width=false,
    variable_dim_scaling=true,
    infill_funcs = [:std,:median]
        )
    
brute = false


# # Optimize the test function
# # Create optimised categorical sampling plan with Categorical(x) possible values in 1:y dimensions
# func = test_funs[:rosenbrock_9D]
# possible_locs = 20
# dims = [Categorical(possible_locs) for i in 1:length(func.sr)]
# vals = [Tuple((i for i in range(-5.0,stop=5.0,length=possible_locs))) for i in 1:length(func.sr)] 

# Optimize the test function
# Create optimised categorical sampling plan with Categorical(x) possible values in 1:y dimensions
func = test_funs[:rosenbrock_2D]
possible_locs = 20
dims = [LatinHypercubeSampling.Categorical(possible_locs) for i in 1:length(func.sr)]
vals = [Tuple((i for i in range(-5.0,stop=5.0,length=possible_locs))) for i in 1:length(func.sr)] 




function closest_index(x_val, vals) 
            
    ibest = first(eachindex(vals)) 
    dxbest = abs(vals[ibest]-x_val) 
    for I in eachindex(vals) 
        dx = abs(vals[I]-x_val) 
        if dx < dxbest 
            dxbest = dx 
            ibest = I 
        end     
    end 
    ibest 
end 

function categorical_smoptimize(func,options,dims,brute,vals)

    #Print the known minimum (bruteforce approach)
    possible_designs = collect(Base.Iterators.product(vals...))    
    println("Number of possible designs = ", length(possible_designs))
    min_loc = argmin(func.fun.(possible_designs))
    max_loc = argmax(func.fun.(possible_designs))
    println("Minimum and maximum value in categorical design space:")
    println("Min ",func.fun(possible_designs[min_loc]), ", ", possible_designs[min_loc])
    println("Max ",func.fun(possible_designs[max_loc]), ", ", possible_designs[max_loc])

    #Load some optional argument values
    @unpack num_start_samples, sampling_plan_opt_gens,
    iterations, trace = options

    search_range = func.sr
    
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
    lhc_samples = SurrogateModelOptim.f_opt_eval(func.fun,lhc_plan,trace)
    criteria = 1

    #Initialize variables to be returned
    sm_interpolant = nothing
    sm_interpolant_cat = nothing
    infill_type = Array{Symbol,1}(undef,0)
    infill_prediction = Array{Float64,1}(undef,0)
    optres = nothing
    infill_plan = Array{Float64,2}(undef,size(lhc_plan,1),0)
    infill_sample = Array{Float64,2}(undef,1,0)

    #Run the optimization iterations number of times
    for i = 1:iterations
        if trace
            print("\n \n \n \t Iteration ")
            printstyled(i,bold=true)
            print(" out of ", iterations, "\n")
        end

        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = SurrogateModelOptim.surrogate_model(samples_all, plan_all; options=options)

        
        sm_interpolant_cat = function (x,vals)
            xc = copy(x)
            for (i,x_val) in enumerate(xc)
                c_ind = closest_index(x_val, vals[i])
                xc[i] = vals[i][c_ind]
            end
            return sm_interpolant(xc)
        end

        #Points to add to the sampling plan to improve the interpolant        
        infill_plan_new, criteria, infill_type_new, infill_prediction_new  = 
        if brute
            SurrogateModelOptim.model_infill_brute(search_range,plan_all,samples_all,x->sm_interpolant_cat(x,vals),vals,criteria;options=options)
        else
            SurrogateModelOptim.model_infill(search_range,plan_all,samples_all,x->sm_interpolant_cat(x,vals),vals,criteria;options=options)
        end

        #Evaluate the new infill points
        infill_sample_new = SurrogateModelOptim.f_opt_eval(func.fun,infill_plan_new,samples_all,trace)

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [infill_sample infill_sample_new]
        infill_type = [infill_type; infill_type_new]
        infill_prediction = [infill_prediction; infill_prediction_new]
    end   

    return SurrogateModelOptim.SurrogateResult( lhc_samples, lhc_plan,                                   #x->sm_interpolant_cat(x,vals),
                    sm_interpolant,
                    optres, infill_sample, infill_type,
                    infill_plan, infill_prediction,options)
end


# This runs num_start_samples + (iterations*num_infill_points) function
# evaluations in total.
result = categorical_smoptimize(func,options,dims,brute,vals)


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
