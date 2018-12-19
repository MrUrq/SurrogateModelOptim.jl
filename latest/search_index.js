var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SurrogateModelOptim.RBFHypersResult",
    "page": "Home",
    "title": "SurrogateModelOptim.RBFHypersResult",
    "category": "type",
    "text": "RBFHypersResult(width::S,kernelFunc,scaling::U,fitness::Float64)\n\nDatastructure to store results from the optimisation of an RBF interpolation kernel\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim.RMSErrorLOO-Tuple{Any,Any,Any}",
    "page": "Home",
    "title": "SurrogateModelOptim.RMSErrorLOO",
    "category": "method",
    "text": "function RMSErrorLOO(interp,samples::Array{Float64,2},plan;\ncond_max = 1/eps(Float64)/10000, rippa = false)\n\nCalculate the Leave-One-Out RMS error for a interpolation method.  rippa can be used to calculate the approximation of the LOO error  for Radial Basis Functions at a cost of 𝛰(3) (compared to 𝛰(4)). cond_max sets the maximum allowed condition number for matrix A used in the RBF calculation.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._rippa-Tuple{Any,Any}",
    "page": "Home",
    "title": "SurrogateModelOptim._rippa",
    "category": "method",
    "text": "function _rippa(A,a)\n\nEstimate the Leave-One-Out (LOO) errors using rippas method. Complexity of 𝛰(3) compared to calculating the exact LOO at a cost of 𝛰(4). A is the RBF matrix and a is the weights of the RBF.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._scale-Union{Tuple{T}, Tuple{Array{T,1},Any,Any}} where T<:Real",
    "page": "Home",
    "title": "SurrogateModelOptim._scale",
    "category": "method",
    "text": "function _scale(oldX::Array{T,1},newMin,newMax;oldMin=minimum(oldX), oldMax=maximum(oldX)) where T <: Real\n\nScale a vector to match a new range.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._scale-Union{Tuple{T}, Tuple{Array{T,2},Int64,Any,Any}} where T<:Real",
    "page": "Home",
    "title": "SurrogateModelOptim._scale",
    "category": "method",
    "text": "function _scale(oldX::Array{T,2},direction::Int,newMin,newMax; oldMin=minimum(oldX,direction)::Array{T,2},oldMax=maximum(oldX,direction)::Array{T,2}) where T <: Real \n\nScale a 2D matrix to match a new range along the specified direction.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._scale-Union{Tuple{T}, Tuple{T,Any,Any}} where T<:Real",
    "page": "Home",
    "title": "SurrogateModelOptim._scale",
    "category": "method",
    "text": "function _scale(oldX::T,newMin,newMax;oldMin=minimum(oldX), oldMax=maximum(oldX)) where T <: Real\n\nScale a scalar to match a new range.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._surrogate_interpolant-Union{Tuple{T}, Tuple{U}, Tuple{T,Any,Any,Any,Any,Any}} where T<:SurrogateModelOptim.RBFHypersResult{U,Float64} where U",
    "page": "Home",
    "title": "SurrogateModelOptim._surrogate_interpolant",
    "category": "method",
    "text": "_surrogate_interpolant(optres::T,points,observations,estimationpoints,\noldMin,oldMax) where T <: RBFoptim_v1.HypersResult{U,Float64} where U\n\nEvaluate an optimised interpolant at locations estimationpoints.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._surrogate_interpolant-Union{Tuple{T}, Tuple{U}, Tuple{T,Any,Any,Any}} where T<:SurrogateModelOptim.RBFHypersResult{U,Float64} where U",
    "page": "Home",
    "title": "SurrogateModelOptim._surrogate_interpolant",
    "category": "method",
    "text": "_surrogate_interpolant(optres::T,plan,samples,estimationpoints,\noldMin,oldMax) where T <: SurrogateModelOptim.RBFHypersResult{U,Float64} where U\n\nEvaluate an optimised interpolant at locations estimationpoints.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim._update_options-Tuple{Any}",
    "page": "Home",
    "title": "SurrogateModelOptim._update_options",
    "category": "method",
    "text": "_update_options(search_range;kwargs...)\n\nInternal function of SurrogateModelOptim to update the options to suitable default values. \n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim.interp_obj-Tuple{Array{Float64,1},Any,Any,Array{Float64,2}}",
    "page": "Home",
    "title": "SurrogateModelOptim.interp_obj",
    "category": "method",
    "text": "function interp_obj(inpt::Vector{Float64}, kerns, samples::Array{Float64,2},\nplan::Array{Float64,2}; rippa::Bool = false, scale::Bool = false)\n\nObjective function for optimisation of interpolation kernel function and width.\n\n\n\n\n\n"
},

{
    "location": "index.html#SurrogateModelOptim.jl-1",
    "page": "Home",
    "title": "SurrogateModelOptim.jl",
    "category": "section",
    "text": "Modules = [SurrogateModelOptim]"
},

]}
