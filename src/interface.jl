function options(;kwargs...)
    kwargs
end

"""
    _update_options(search_range;kwargs...)

Internal function of SurrogateModelOptim to update the options to suitable default
values. 
"""
function _update_options(search_range;kwargs...)
    
    #Instantiate default options list
    def_options = Options()     
    
    #Update options with user settings
    options = reconstruct(def_options,kwargs)                   
        
    #Set the number of sample points as 10 times the number of dimensions by default
    if :num_start_samples ∉ keys(kwargs)
        options = reconstruct(options,num_start_samples = length(search_range)*10) 
    end    

    #No dimensional scaling for 1-D problems unless specified by the user
    if (:variable_dim_scaling ∉ keys(kwargs)) && (length(search_range) == 1)
        options = reconstruct(options,variable_dim_scaling = false)  
    end 

    #Make sure the correct options for smoothing is used
    @assert ((options.smooth == false) || (options.smooth == :variable) || 
     (options.smooth == :single) || (options.smooth == :single_user)) "Not supported option for smooth"

    return options
end


function _LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens)
    unscaled_plan = LHCoptim(num_start_samples,length(search_range),sampling_plan_opt_gens)[1]'

    plan = Array{eltype(eltype(search_range))}(undef,size(unscaled_plan))
    for (i,sr) in enumerate(search_range)
        plan[i,:] = _scale(unscaled_plan[i,:],sr[1],sr[2])
    end
    return plan
end


            