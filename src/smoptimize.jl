"""
    smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}; options=Options())

Optimize the function `f` in the range `search_range` using a Radial Basis Function based surrogate model.
"""
function smoptimize(f::Function, search_range::Array{Tuple{Float64,Float64},1}; options::Options=Options())

    #Load some option values
    @unpack num_start_samples, sampling_plan_opt_gens,
            iterations, trace, create_final_surrogate = options
    
    #Create sampling plan
    lhc_plan = scaled_LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens;trace=trace)
    
    #Evaluate sampling plan
    lhc_samples = f_opt_eval(f,lhc_plan;trace=trace)

    #Initialize variables to be returned
    sm_interpolant = nothing
    infill_type = Array{Symbol,1}(undef,0)
    infill_prediction = Array{Float64,1}(undef,0)
    optres = nothing
    infill_plan = Array{Float64,2}(undef,size(lhc_plan,1),0)
    infill_sample = Array{Float64,2}(undef,1,0)

    #Run the entire optimization iterations number of times
    for i = 1:iterations
        
        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = surrogate_model(plan_all, samples_all; options=options)
        
        #Points to add to the sampling plan to improve the interpolant
        infill_plan_new, infill_type_new, infill_prediction_new, options  = model_infill(search_range,plan_all,
                                                                                samples_all,sm_interpolant;options=options)
        
        #Evaluate the new infill points
        print_iteration(trace,i,iterations)
        infill_sample_new = f_opt_eval(f,infill_plan_new,samples_all;trace=trace)
        

        #Add infill points
        infill_plan = [infill_plan infill_plan_new]
        infill_sample = [infill_sample infill_sample_new]
        infill_type = [infill_type; infill_type_new]
        infill_prediction = [infill_prediction; infill_prediction_new]

    end

    if create_final_surrogate
        #Create the optimized Radial Basis Function interpolant      
        samples_all = [lhc_samples infill_sample]
        plan_all = [lhc_plan infill_plan]
        sm_interpolant, optres = surrogate_model(plan_all, samples_all; options=options)
    end
    
    return SurrogateResult( lhc_samples, lhc_plan, sm_interpolant,
                            optres, infill_sample, infill_type,
                            infill_plan, infill_prediction,options)
end