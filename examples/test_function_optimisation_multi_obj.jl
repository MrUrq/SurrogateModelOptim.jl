using PlotlyJS
using Statistics
using SurrogateModelOptim
using BlackBoxOptim
using Parameters
using Printf

options=SurrogateModelOptim.Options(
    iterations=10,
    smooth=false, # Should only be false for smooth deterministic functions (noise free)
    num_interpolants=10, #Preferably even number of added processes
    num_start_samples=4,
    constrained_seed_gens=0,
    rbf_opt_gens=1000,
    create_final_surrogate=true, #Use the results from last iteration to
                                 #re-create the surrogate before using it for plotting
        )

fitness_schaffer1(x) = (sum(abs2, x), sum(abs2, x .- 2.0))

search_range = [(-10.0,10.0) for _ in 1:3]


function multi_objective_smoptimize(f,search_range,options)

    #2 objectives
    n_funcs = 2
    
    #Load some option values
    @unpack num_start_samples, sampling_plan_opt_gens,
            iterations, trace = options

    lhc_plan = SurrogateModelOptim.scaled_LHC_sampling_plan(search_range,num_start_samples,10000)
       
    #Evaluate sampling plan
    lhc_samples = [SurrogateModelOptim.f_opt_eval(x->f(x)[i],lhc_plan) for i in 1:n_funcs]

    #Initialize variables to be returned
    sm_interpolant = [nothing for _ in 1:n_funcs]
    infill_type = [Array{Symbol,1}(undef,0) for _ in 1:n_funcs]
    infill_prediction = [Array{Float64,1}(undef,0) for _ in 1:n_funcs]
    optres = [nothing for _ in 1:n_funcs]
    infill_plan = Array{Float64,2}(undef,size(lhc_plan,1),0)
    infill_sample = [Array{Float64,2}(undef,1,0) for _ in 1:n_funcs]
    sm_med = nothing
    

    #Run the entire optimization iterations number of times
    for i = 1:iterations
        SurrogateModelOptim.print_iteration(trace,i,iterations)
        
        plan_all = [lhc_plan infill_plan]
        samples_all = [[lhc_samples[j] infill_sample[j]] for j in 1:n_funcs]
        
        #Create the optimized Radial Basis Function interpolant for each objective function
        sm_fs = Tuple(SurrogateModelOptim.surrogate_model(plan_all, samples_all[j]; options=options)[1] for j in 1:n_funcs)
        sm = x -> Tuple(sm_fs[i](x) for i in 1:n_funcs)
    
        #We add two points per iteration, one exploring and one exploiting sample
        sm_med = x-> median.(sm(x))
        sm_std = x-> -1 .*std.(sm(x))
        
        res_med = bboptimize(sm_med; Method=:borg_moea,
                FitnessScheme=ParetoFitnessScheme{n_funcs}(is_minimizing=true),
                SearchRange=search_range, NumDimensions=size(plan_all,2), ϵ=0.05,
                MaxSteps=10_000, TraceInterval=1.0, TraceMode=:silent);

        res_std = bboptimize(sm_std; Method=:borg_moea,
                FitnessScheme=ParetoFitnessScheme{n_funcs}(is_minimizing=true),
                SearchRange=search_range, NumDimensions=size(plan_all,2), ϵ=0.05,
                MaxSteps=10_000, TraceInterval=1.0, TraceMode=:silent);
        

        #Add two points randomly selected points from the pareto front
        pf_res_med = pareto_frontier(res_med); p_med = pf_res_med[rand(1:length(pf_res_med))]
        pf_res_std = pareto_frontier(res_std); p_std = pf_res_std[rand(1:length(pf_res_std))]

        infill_plan_new = [permutedims(p_med.inner.params') permutedims(p_std.inner.params')]
            
        #Evaluate the new infill points
        infill_sample_new = [SurrogateModelOptim.f_opt_eval(x->f(x)[j],infill_plan_new,samples_all[j];trace=:silent) for j in 1:n_funcs]

        println("\nEvaluating functions------------")
        [str = printstyled(@sprintf("Objective value: %-15.7g \n",minimum(infill_sample_new[i]))) for i in 1:n_funcs]
        println("--------------------------------")

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [[infill_sample[j] infill_sample_new[j]] for j in 1:n_funcs]
        samples_all = [[lhc_samples[j] infill_sample[j]] for j in 1:n_funcs]
    end   
    
    return sm_med
end
sm = multi_objective_smoptimize(fitness_schaffer1,search_range,options)




function plot_pareto(fun,sm)   
    
    res_bbopt = bboptimize(fun; Method=:borg_moea,
                FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                SearchRange=search_range, NumDimensions=3, ϵ=0.05,
                MaxSteps=50000, TraceInterval=1.0, TraceMode=:silent);
    pf1= pareto_frontier(res_bbopt)
    bb_f1 = [pf1[i].fitness[1] for i in 1:length(pf1)]
    bb_f2 = [pf1[i].fitness[2] for i in 1:length(pf1)]

    res_smopt = bboptimize(sm; Method=:borg_moea,
                FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                SearchRange=search_range, NumDimensions=3, ϵ=0.05,
                MaxSteps=50000, TraceInterval=1.0, TraceMode=:silent);
    pf2 = pareto_frontier(res_smopt)
    sm_f1 = [pf2[i].fitness[1] for i in 1:length(pf2)]
    sm_f2 = [pf2[i].fitness[2] for i in 1:length(pf2)]

    trace1 = scatter(;x=bb_f1, y=bb_f2, mode="markers",name= "BB_optim")
    trace2 = scatter(;x=sm_f1, y=sm_f2, mode="markers",name= "SM_optim")
    layout = Layout(title="Pareto")
    plot([trace1, trace2],layout)
end

# Plot the results 
display(plot_pareto(fitness_schaffer1,sm))