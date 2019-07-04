"""
    scaled_LHC_sampling_plan(search_range::Array{Tuple{Float64,Float64},1},num_samples,sampling_plan_opt_gens;trace::Symbol=:silent)

Create a Latin Hypercube sampling plan scaled to the `search_range` where 
`size(plan) = (num_dimensions,num_samples)`.
"""
function scaled_LHC_sampling_plan(search_range::Array{Tuple{Float64,Float64},1},num_samples,sampling_plan_opt_gens;trace::Symbol=:silent)
   
    #Printing of creation
    (trace == :verbose) && print_lhc_create()

    unscaled_plan, audze = LHCoptim(num_samples,length(search_range),sampling_plan_opt_gens)

    unscaled_plan = unscaled_plan'

    #Verbose printing of resulting plan
    (trace == :verbose) && print_lhc_plan(audze[end],num_samples,length(search_range),sampling_plan_opt_gens)

    plan = Array{eltype(eltype(search_range))}(undef,size(unscaled_plan))
    for (i,sr) in enumerate(search_range)
        plan[i,:] = scale(unscaled_plan[i,:],sr[1],sr[2])
    end
    return plan
end          

function print_lhc_create()
    println("Creating optimized Latin Hypercube Sampling Plan")    
end

function print_lhc_plan(audze,num_samples,dims,sampling_plan_opt_gens)
    println("    Optimised LHC plan summary")
    println(@sprintf("    %-15s %-40s %-10s %-0s", "Iterations","Audze–Eglais obj. (higher is better)","Samples","Dimensions"))
    println(@sprintf("    %-15.7g %-40.7g %-10.7g %-0.7g \n", sampling_plan_opt_gens,audze,num_samples,dims))
end