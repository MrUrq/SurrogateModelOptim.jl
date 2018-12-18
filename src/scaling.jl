"""
function _scale(oldX::T,newMin,newMax;oldMin=minimum(oldX),
oldMax=maximum(oldX)) where T <: Real

Scale a scalar to match a new range.
"""
function _scale(oldX::T,newMin,newMax;oldMin=minimum(oldX),oldMax=maximum(oldX)) where T <: Real

oldX = (((oldX - oldMin) * (newMax - newMin)) / (oldMax - oldMin)) + newMin
end

"""
function _scale(oldX::Array{T,1},newMin,newMax;oldMin=minimum(oldX),
oldMax=maximum(oldX)) where T <: Real

Scale a vector to match a new range.
"""
function _scale(oldX::Array{T,1},newMin,newMax;oldMin=minimum(oldX),oldMax=maximum(oldX)) where T <: Real

newX = (((oldX .- oldMin) .* (newMax .- newMin)) ./ (oldMax .- oldMin)) .+ newMin
end

"""
function _scale(oldX::Array{T,2},direction::Int,newMin,newMax;
oldMin=minimum(oldX,direction)::Array{T,2},oldMax=maximum(oldX,direction)::Array{T,2}) where T <: Real 

Scale a 2D matrix to match a new range along the specified `direction`.
"""
function _scale(oldX::Array{T,2},direction::Int,newMin,newMax;oldMin=minimum(oldX,dims=direction)::Array{T,2},oldMax=maximum(oldX,dims=direction)::Array{T,2}) where T <: Real 

newX = mapslices(x -> _scale(x,newMin,newMax), oldX, dims=direction)
end