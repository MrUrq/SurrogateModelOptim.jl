"""
    scaled_LHC_sampling_plan(search_range::Array{Tuple{Float64,Float64},1},num_samples,sampling_plan_opt_gens;trace=false)

Create a Latin Hypercube sampling plan scaled to the `search_range` where 
`size(plan) = (num_dimensions,num_samples)`.
"""
function scaled_LHC_sampling_plan(search_range::Array{Tuple{Float64,Float64},1},num_samples,sampling_plan_opt_gens;trace=false)
   
    #Printing of creation
    trace && println("Creating optimized Latin Hypercube Sampling Plan ... \n")

    unscaled_plan = LHCoptim(num_samples,length(search_range),sampling_plan_opt_gens)[1]'

    plan = Array{eltype(eltype(search_range))}(undef,size(unscaled_plan))
    for (i,sr) in enumerate(search_range)
        plan[i,:] = _scale(unscaled_plan[i,:],sr[1],sr[2])
    end
    return plan
end          