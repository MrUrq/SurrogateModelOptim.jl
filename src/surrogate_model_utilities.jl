"""
    function _rippa(A,a)

Estimate the Leave-One-Out (LOO) errors using rippas method. Complexity of 𝛰(3)
compared to calculating the exact LOO at a cost of 𝛰(4). `A` is the RBF matrix
and `a` is the weights of the RBF.
"""
function _rippa(A, a)
    N = size(A, 1)           #Number of Leave-One-Out (LOO) errors
    e = Matrix(I, N, N)     
    E = Array{Float64}(undef,N)

    for k = 1:N
        xₖ = A \ e[:,k]       #xₖ solution to Ax[k] = e[k]   - Solved N times
        E[k] = a[k] / xₖ[k]   #Estimated error for k-th subset
    end

    return E
end

"""
    function RMSErrorLOO(interp,samples::Array{Float64,2},plan;
    cond_max = 1/eps(Float64)/10000, rippa = false)

Calculate the Leave-One-Out RMS error for a interpolation method. 
`rippa` can be used to calculate the approximation of the LOO error 
for Radial Basis Functions at a cost of 𝛰(3)
(compared to 𝛰(4)). `cond_max` sets the maximum allowed condition number for
matrix `A` used in the RBF calculation.
"""
function RMSErrorLOO(interp, samples, plan, smooth;
    cond_max::Float64=1e6, rippa::Bool=false, rbf_dist_metric = Euclidean())

    #initiate arrays
    N = length(samples)
    E = Array{Float64}(undef,N)       #Error
    A = zeros(Float64, N, N)    #RBF matrix A
    ests = Array{Float64}(undef,N)    #Function estimate based on LOO
    LOOinds = Array{Int}(undef,N - 1, N)
    for i = 1:N
        LOOinds[:,i] = filter(x -> x != i, 1:N) # Get the leave one out sub indices
    end

    RMSErrorLOO!(E, A, ests, LOOinds, interp, samples, plan, smooth;
    cond_max=cond_max, rippa=rippa, rbf_dist_metric = rbf_dist_metric)
end


function RMSErrorLOO!(E, A, ests, LOOinds, interp::U,
 samples, plan::T, smooth; cond_max::Float64=1e6,
  rippa::Bool=false, rbf_dist_metric = Euclidean()) where T <: AbstractArray where U <: AbstractArray

    N = length(samples)

    #perform Leave-One-Out estimation and calculate error
    if rippa
        @assert typeof(interp) <: Vector{ScatteredInterpolation.RadialBasisFunction} "Rippas
              algorithm only available for Radial Basis Functions"

        #interpolation object trained on the entire dataset
        itp, A = interpolate(interp, plan,
            samples, returnRBFmatrix=true, metric = rbf_dist_metric, smooth=smooth)

        #RBF error estimation based on Rippas algorithm
        E = _rippa(A, itp.w)

    else
        for i = 1:N

            if typeof(smooth) <: AbstractVector
                loo_smooth = smooth[LOOinds[:,i]]
            else
                loo_smooth = smooth
            end

            #interpolation object trained on the LOO information
            itp, A = interpolate(interp[LOOinds[:,i]], plan[:,LOOinds[:,i]],
                samples[LOOinds[:,i]], returnRBFmatrix=true, metric = rbf_dist_metric, smooth=loo_smooth)

             #evaluate the interpolation object in the LOO position
            ests[i] = ScatteredInterpolation.evaluate(itp, plan[:,i])[1]

            E[i] = samples[i] - ests[i]
        end
    end

    #check the conditioning of matrix A for RBFs
    isrbf = typeof(interp) <: ScatteredInterpolation.RadialBasisFunction

    if isrbf && (cond(A) > cond_max)
        RMSE = Inf
    else
        RMSE = sqrt(mean(E.^2))
    end
    return RMSE
end




function RMSErrorLOO!(E, A, ests, LOOinds, interp, samples, plan::T, smooth;
 cond_max::Float64=1e6, rippa::Bool=false, rbf_dist_metric = Euclidean()) where T <: AbstractArray

    N = length(samples)

    #perform Leave-One-Out estimation and calculate error
    if rippa
        @assert typeof(interp) <: ScatteredInterpolation.RadialBasisFunction "Rippas algorithm only available for Radial Basis Functions"

        #interpolation object trained on the entire dataset
        itp, A = interpolate(interp, plan,
         samples, returnRBFmatrix=true, metric = rbf_dist_metric, smooth=smooth)

        #RBF error estimation based on Rippas algorithm
        E = _rippa(A, itp.w)

    else
        for i = 1:N

            if typeof(smooth) <: AbstractVector
                loo_smooth = smooth[LOOinds[:,i]]
            else
                loo_smooth = smooth
            end

            #interpolation object trained on the LOO information
            itp, A = interpolate(interp, plan[:,LOOinds[:,i]],
                samples[LOOinds[:,i]], returnRBFmatrix=true, metric = rbf_dist_metric, smooth=loo_smooth)

            #evaluate the interpolation object in the LOO position
            ests[i] = ScatteredInterpolation.evaluate(itp, plan[:,i])[1]

            E[i] = samples[i] - ests[i]
        end
    end

    #check the conditioning of matrix A for RBFs
    isrbf = typeof(interp) <: ScatteredInterpolation.RadialBasisFunction

    if isrbf && (cond(A) > cond_max)
        RMSE = Inf
    else
        RMSE = sqrt(mean(E.^2))
    end
    return RMSE
end

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

_scale(x::Missing,min_val,max_val;old_min,old_max) = missing


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

function preprocess_point!(preprocessed_point,points,optres;old_min::Array{Float64,2},old_max::Array{Float64,2})

    for i = 1:size(preprocessed_point,1)
        preprocessed_point[i,:] = _scale(points[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i],
        old_min = old_min[i], old_max = old_max[i])
    end
    return preprocessed_point
end

function preprocess_point!(preprocessed_point,points,optres::SurrogateModelOptim.RBFHypers{Bool,U};old_min::Array{Float64,2},old_max::Array{Float64,2}) where U
    preprocessed_point = points
end

function preprocess_point(points,optres::SurrogateModelOptim.RBFHypers{Bool,U};base_scale::Array{Float64,2}) where U
    points    
end

function construct_search_range(plan::Array{Float64,2}, variable_kernel_width,
    min_rbf_width, max_rbf_width, variable_dim_scaling,
    min_scale, max_scale, smooth, max_smooth)

n_dims, n_samples = size(plan)
sr = Array{Tuple{Float64,Float64},1}()


# Add the kernel width and type search range
variable_kernel_width ? n_kerns = n_samples : n_kerns = 1    
push!(sr,   create_sr(  (min_range=min_rbf_width, max_range=max_rbf_width, n_times=n_kerns),
(min_range=0.0, max_range=1.0, n_times=n_kerns))...)

# Add the dimensional scaling search range
variable_dim_scaling ? n_dim_scales = n_dims : n_dim_scales = 0
push!(sr,   create_sr(  (min_range=min_scale, max_range=max_scale, n_times=n_dim_scales))...)

# Add the ridge regression smoothing search range
(smooth == :variable) && (n_smooth = n_samples)
(smooth == :single)   && (n_smooth = 1)
(smooth == :single_user)   && (n_smooth = 0)
(smooth == false)     && (n_smooth = 0)
push!(sr,   create_sr(  (min_range=0.0, max_range=max_smooth, n_times=n_smooth))...)

return sr
end

function create_sr(vargs::NamedTuple{(:min_range, :max_range, :n_times),Tuple{Float64,Float64,Int64}}...)
sr = Array{Tuple{Float64,Float64},1}()
for varg in vargs
for i = 1:varg[3]
push!(sr,(varg[1],varg[2]))
end
end        
return sr
end

function extract_vector_range(vargs::Int64...)    
output = Array{Any,1}()

count = 0
for (i,varg) in enumerate(vargs)
if varg == 0
push!(output, false)
else
push!(output, (count+1):(count+varg))
count += varg      
(length(output[end]) == 1) && (output[end] = output[end][1])
end
end
return Tuple(output)
end

function extract_bboptim_hypers(bboptim_fcall_vector,plan,kerns,
    variable_kernel_width,variable_dim_scaling,
    smooth,smooth_user)

n_dims, n_samples = size(plan)

# Kernel width and type length
variable_kernel_width ? n_kerns = n_samples : n_kerns = 1    

# Dimensional scaling length
variable_dim_scaling ? n_dim_scales = n_dims : n_dim_scales = 0

# Ridge regression smoothing length
(smooth == :variable) && (n_smooth = n_samples)
(smooth == :single)   && (n_smooth = 1)
(smooth == :single_user)   && (n_smooth = 0)
(smooth == false)     && (n_smooth = 0)

width_inds, kernel_float_inds, scaling_inds, smooth_inds = extract_vector_range(n_kerns,
                                                        n_kerns,
                                                        n_dim_scales,
                                                        n_smooth
                                                        )

# Arrange the width and smoothing results
width = bboptim_fcall_vector[width_inds]
kernel_float = bboptim_fcall_vector[kernel_float_inds]
!(scaling_inds == false) ? scaling = bboptim_fcall_vector[scaling_inds] : scaling = scaling_inds
(smooth == false) &&            (smooth = smooth_inds)
(smooth == :single) &&          (smooth = bboptim_fcall_vector[smooth_inds])
(smooth == :variable) &&        (smooth = bboptim_fcall_vector[smooth_inds])
(smooth == :single_user) &&     (smooth = smooth_user)



# Arrange the RBF kernel result
kern_ind = round.(Int,_scale(kernel_float,1,length(kerns),old_min=0,old_max=1))
if variable_kernel_width
kern = Vector{ScatteredInterpolation.RadialBasisFunction}(undef,size(plan,2))
for i = 1:size(plan,2)
kern[i] = kerns[kern_ind[i]](width[i])
end
elseif !variable_kernel_width
kern = kerns[kern_ind](width)
end

return kern, scaling, smooth
end

