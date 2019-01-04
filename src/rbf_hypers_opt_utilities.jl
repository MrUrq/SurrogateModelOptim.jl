function _interp_fixK_fixW_obj(inpt::Vector{Float64}, kerns, samples,
    plan::Array{Float64,2}; rippa::Bool = false, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)

    x = inpt[1]
    y = inpt[2]

    # Choose interpolation kernel function based on floating point
    kern_ind = round.(Int,_scale(y,1,length(kerns),old_min=0,old_max=1))
    kern = kerns[kern_ind](x)

    try
        E = RMSErrorLOO(kern, samples, plan; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
    catch
        E = Inf
    end    
end


function _scaled_interp_fixK_fixW_obj(inpt::Vector{Float64}, kerns, samples,
    plan::Array{Float64,2}; rippa::Bool = false, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)

    x = inpt[1]     #Kernel width
    y = inpt[2]     #Kernel function
    z = inpt[3:end] #Axis scaling

    @assert length(z) == size(plan,1)   #make sure all dimensions are scaled

    # Scale the plan 
    scaled_plan = similar(plan)
    for i = 1:size(plan,1)
        scaled_plan[i,:] = _scale(plan[i,:],-1.0*z[i],1.0*z[i])
    end

    # Choose interpolation kernel function based on floating point
    kern_ind = round.(Int,_scale(y,1.0,float(length(kerns)),old_min=0,old_max=1))
    kern = kerns[kern_ind](x)

    try
        E = RMSErrorLOO(kern, samples, scaled_plan; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
    catch
        E = Inf
    end
end


function _interp_varK_varW_obj(inpt::Vector{Float64}, kerns, samples,
    plan::Array{Float64,2}; rippa::Bool = false, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)

    x = inpt[1:size(plan,2)]     #Kernel width
    y = inpt[size(plan,2)+1:size(plan,2)+size(plan,2)] #Kernel function

    # Choose interpolation kernel function based on floating point
    kern_ind = round.(Int,_scale(y,1.0,float(length(kerns)),old_min=0,old_max=1))
    kern = Vector{ScatteredInterpolation.RadialBasisFunction}(undef,size(plan,2))
    for i = 1:size(plan,2)
        kern[i] = kerns[kern_ind[i]](x[i])
    end

    try
        E = RMSErrorLOO(kern, samples, plan; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
    catch 
        E = Inf
    end
end


function _scaled_interp_varK_varW_obj(inpt::Vector{Float64}, kerns, samples,
    plan::Array{Float64,2}; rippa::Bool = false, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)

    x = inpt[1:size(plan,2)]     #Kernel width
    y = inpt[size(plan,2)+1:size(plan,2)+size(plan,2)] #Kernel function
    z = inpt[size(plan,2)+size(plan,2)+1:end] #Axis scaling


    @assert length(z) == size(plan,1)   #make sure all dimensions are scaled

    # Scale the plan 
    scaled_plan = similar(plan)
    for i = 1:size(plan,1)
        scaled_plan[i,:] = _scale(plan[i,:],-1.0*z[i],1.0*z[i])
    end
    

    # Choose interpolation kernel function based on floating point
    kern_ind = round.(Int,_scale(y,1.0,float(length(kerns)),old_min=0,old_max=1))
    kern = Vector{ScatteredInterpolation.RadialBasisFunction}(undef,size(plan,2))
    for i = 1:size(plan,2)
        kern[i] = kerns[kern_ind[i]](x[i])
    end

    try
        E = RMSErrorLOO(kern, samples, scaled_plan; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
    catch
        E = Inf
    end
end


function point_scale(points,optres;base_scale::Array{Float64,2})

    old_min = minimum(base_scale,dims=2)
    old_max = maximum(base_scale,dims=2)

    scaled_points = similar(points)
    for i = 1:size(scaled_points,1)
        scaled_points[i,:] = _scale(points[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i],
        old_min = old_min[i], old_max = old_max[i])
    end
    return scaled_points
end



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
function RMSErrorLOO(interp, samples, plan;
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

    RMSErrorLOO!(E, A, ests, LOOinds, interp, samples, plan;
    cond_max=cond_max, rippa=rippa, rbf_dist_metric = rbf_dist_metric)
end


function RMSErrorLOO!(E, A, ests, LOOinds, interp::U,
 samples, plan::T; cond_max::Float64=1e6,
  rippa::Bool=false, rbf_dist_metric = Euclidean()) where T <: AbstractArray where U <: AbstractArray

    N = length(samples)

    #perform Leave-One-Out estimation and calculate error
    if rippa
        @assert typeof(interp) <: Vector{ScatteredInterpolation.RadialBasisFunction} "Rippas
              algorithm only available for Radial Basis Functions"

     #interpolation object trained on the entire dataset
        itp, A = interpolate(interp, plan,
      samples, returnRBFmatrix=true, metric = rbf_dist_metric)

     #RBF error estimation based on Rippas algorithm
        E = _rippa(A, itp.w)

    else
        for i = 1:N
         #interpolation object trained on the LOO information
            itp, A = interpolate(interp[LOOinds[:,i]], plan[:,LOOinds[:,i]],
             samples[LOOinds[:,i]], returnRBFmatrix=true, metric = rbf_dist_metric)

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




function RMSErrorLOO!(E, A, ests, LOOinds, interp, samples, plan::T;
 cond_max::Float64=1e6, rippa::Bool=false, rbf_dist_metric = Euclidean()) where T <: AbstractArray

    N = length(samples)

    #perform Leave-One-Out estimation and calculate error
    if rippa
        @assert typeof(interp) <: ScatteredInterpolation.RadialBasisFunction "Rippas algorithm only available for Radial Basis Functions"

        #interpolation object trained on the entire dataset
        itp, A = interpolate(interp, plan,
         samples, returnRBFmatrix=true, metric = rbf_dist_metric)

        #RBF error estimation based on Rippas algorithm
        E = _rippa(A, itp.w)

    else
        for i = 1:N
        #interpolation object trained on the LOO information
        itp, A = interpolate(interp, plan[:,LOOinds[:,i]],
         samples[LOOinds[:,i]], returnRBFmatrix=true, metric = rbf_dist_metric)

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
