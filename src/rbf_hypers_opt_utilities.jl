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



