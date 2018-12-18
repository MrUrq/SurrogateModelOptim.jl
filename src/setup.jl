function _interface_checks(SearchRange, NumDimensions)
    #Correct SearchRange type signature
    @assert eltype(SearchRange) == Float64 || 
    eltype(SearchRange) == Tuple{Float64,Float64} "Using $(typeof(SearchRange)) for SearchRange is not supported."

    #Correct use of the dimensionality and search range.
    @assert ((eltype(SearchRange) == Float64) ? (typeof(NumDimensions) == Int64) : true) "Specify the dimensionality of the problem with keyword argument NumDimensions"
    @assert ((eltype(SearchRange) == Tuple{Float64,Float64}) ? (typeof(NumDimensions) == Bool) : true) "Cannot specify NumDimensions and provide SearchRange as an array of ranges"
end

"""
    _setup_parameters(SearchRange, NumDimensions=false; kwargs...)

Internal function of SurrogateModelOptim to handle the interface which constrains 
the input combinations the user can use. 
"""
function _setup_parameters(SearchRange, NumDimensions=false; kwargs...)

    _interface_checks(SearchRange, NumDimensions)
      
    #Instantiate default parameters list
    params = _Parameters()     
    
    #Update parameters with user settings
    params = Parameters.reconstruct(params,kwargs)                   
        
    #Set the number of sample points as 10 times the number of dimensions by default
    if :NumStartSamples ∉ keys(kwargs)
        if NumDimensions == false
            params = Parameters.reconstruct(params,NumStartSamples = length(SearchRange)*10)  
        else
            params = Parameters.reconstruct(params,NumStartSamples = NumDimensions*10) 
        end
    end    

    return params
end


function _LHC_sampling_plan(SearchRange,NumDimensions::Bool,NumStartSamples,SamplingPlanOptGens)
    unscaled_plan = LatinHypercubeSampling.LHCoptim(NumStartSamples,length(SearchRange),SamplingPlanOptGens)[1]'

    plan = Array{eltype(eltype(SearchRange))}(undef,size(unscaled_plan))
    for (i,sr) in enumerate(SearchRange)
        plan[i,:] = _scale(unscaled_plan[i,:],sr[1],sr[2])
    end
    return plan
end


function _LHC_sampling_plan(SearchRange,NumDimensions::Int,NumStartSamples,SamplingPlanOptGens)
    unscaled_plan = LatinHypercubeSampling.LHCoptim(NumStartSamples,NumDimensions,SamplingPlanOptGens)[1]'

    plan = Array{eltype(eltype(SearchRange))}(undef,size(unscaled_plan))
    for i = 1:NumDimensions
        plan[i,:] = _scale(unscaled_plan[i,:],SearchRange[1],SearchRange[2])
    end
    return plan
end