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