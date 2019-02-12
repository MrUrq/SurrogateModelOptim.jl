
function _LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens,trace)
   
    #Printing of creation
    trace && _LHC_trace()

    unscaled_plan = LHCoptim(num_start_samples,length(search_range),sampling_plan_opt_gens)[1]'

    plan = Array{eltype(eltype(search_range))}(undef,size(unscaled_plan))
    for (i,sr) in enumerate(search_range)
        plan[i,:] = _scale(unscaled_plan[i,:],sr[1],sr[2])
    end
    return plan
end

function _LHC_trace()
    println("Creating optimized Latin Hypercube Sampling Plan ... \n")
end


            