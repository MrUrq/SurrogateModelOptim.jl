"""
    function _scale(old_x::T,new_min,new_max;old_min=minimum(old_x),
    old_max=maximum(old_x)) where T <: Real

Scale a scalar to match a new range.
"""
function _scale(old_x::T,new_min,new_max;old_min=minimum(old_x),old_max=maximum(old_x)) where T <: Real

    old_x = (((old_x - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
end

"""
    function _scale(old_x::Array{T,1},new_min,new_max;old_min=minimum(old_x),
    old_max=maximum(old_x)) where T <: Real

Scale a vector to match a new range.
"""
function _scale(old_x::Array{T,1},new_min,new_max;old_min=minimum(old_x),old_max=maximum(old_x)) where T <: Real

    newX = (((old_x .- old_min) .* (new_max .- new_min)) ./ (old_max .- old_min)) .+ new_min
end

"""
    function _scale(old_x::Array{T,2},direction::Int,new_min,new_max;
    old_min=minimum(old_x,direction)::Array{T,2},old_max=maximum(old_x,direction)::Array{T,2}) where T <: Real 

Scale a 2D matrix to match a new range along the specified `direction`.
"""
function _scale(old_x::Array{T,2},direction::Int,new_min,new_max;old_min=minimum(old_x,dims=direction)::Array{T,2},old_max=maximum(old_x,dims=direction)::Array{T,2}) where T <: Real 

    newX = mapslices(x -> _scale(x,new_min,new_max), old_x, dims=direction)
end


function preprocess_point(points,optres;base_scale::Array{Float64,2})

    old_min = minimum(base_scale,dims=2)
    old_max = maximum(base_scale,dims=2)

    preprocessed_point = similar(points)
    for i = 1:size(preprocessed_point,1)
        preprocessed_point[i,:] = _scale(points[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i],
        old_min = old_min[i], old_max = old_max[i])
    end
    return preprocessed_point
end

function preprocess_point(points,optres::SurrogateModelOptim.RBFHypers{Bool,U};base_scale::Array{Float64,2}) where U
    points    
end