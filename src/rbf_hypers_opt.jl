"""
    function interp_obj(inpt::Vector{Float64}, kerns, samples::Array{Float64,2},
    plan::Array{Float64,2}; rippa::Bool = false, scale::Bool = false)

Objective function for optimisation of interpolation kernel function and width.
"""
function interp_obj(inpt::Vector{Float64}, kerns, samples,
        plan::Array{Float64,2}; rippa::Bool = false,
        variable_kernel_width::Bool = true, variable_dim_scaling::Bool = true,
        cond_max=cond_max, rbf_dist_metric = Distances.Euclidean())

    if variable_kernel_width
        if variable_dim_scaling
            return _scaled_interp_varK_varW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        else
            return _interp_varK_varW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        end
    elseif !variable_kernel_width
        if variable_dim_scaling
            return _scaled_interp_fixK_fixW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        else
            return _interp_fixK_fixW_obj(inpt,kerns,samples,plan;
            rippa = rippa, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
        end        
    end
end


function _interp_fixK_fixW_obj(inpt::Vector{Float64}, kerns, samples,
    plan::Array{Float64,2}; rippa::Bool = false, cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)

    x = inpt[1]
    y = inpt[2]

    # Choose interpolation kernel function based on floating point
    kernind = round.(Int,_scale(y,1,length(kerns),oldMin=0,oldMax=1))
    kern = kerns[kernind](x)

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
    kernind = round.(Int,_scale(y,1.0,float(length(kerns)),oldMin=0,oldMax=1))
    kern = kerns[kernind](x)

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
    kernind = round.(Int,_scale(y,1.0,float(length(kerns)),oldMin=0,oldMax=1))
    kern = Vector{ScatteredInterpolation.RadialBasisFunction}(undef,size(plan,2))
    for i = 1:size(plan,2)
        kern[i] = kerns[kernind[i]](x[i])
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
    kernind = round.(Int,_scale(y,1.0,float(length(kerns)),oldMin=0,oldMax=1))
    kern = Vector{ScatteredInterpolation.RadialBasisFunction}(undef,size(plan,2))
    for i = 1:size(plan,2)
        kern[i] = kerns[kernind[i]](x[i])
    end

    try
        E = RMSErrorLOO(kern, samples, scaled_plan; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric)
    catch
        E = Inf
    end
end



function _rbf_hypers_opt(samples::Array{Float64,2}, plan::Array{Float64,2}, options::SurrogateModelOptim.Options)
    
    @unpack rippa, variable_kernel_width, variable_dim_scaling, rbf_opt_method, 
            max_rbf_width, max_scale, cond_max, rbf_dist_metric,
            rbf_opt_gens, kerns = options

    nFuncs, nObs = size(samples)
    nDims = size(plan,1)

    # optimise hypers using the same kernel and width for each point
    if variable_kernel_width

        # initiate exploration space for optimisation
        sr = [repeat([(1e-4, max_rbf_width)],nObs); repeat([(0.0, 1.0)],nObs)]
        if variable_dim_scaling
            sr = [sr; repeat([(1e-4, max_scale)],nDims)]
            optHypRes = Array{RBFHypersResult{Array{Float64,1},Array{Float64,1}},1}(undef,nFuncs)
        else
            optHypRes = Array{RBFHypersResult{Array{Float64,1},Float64},1}(undef,nFuncs)            
        end    

    # optimise hypers using different kernels and different widths
    elseif !variable_kernel_width
                
        # initiate exploration space for optimisation
        sr = [(1e-4, max_rbf_width), (0.0, 1.0)]
        if variable_dim_scaling 
            sr = [sr; repeat([(1e-4, max_scale)],nDims)]
            optHypRes = Array{RBFHypersResult{Float64,Array{Float64,1}},1}(undef,nFuncs)
        else
            optHypRes = Array{RBFHypersResult{Float64,Float64},1}(undef,nFuncs)
        end
    
    else
        error("not supported combination of inputs")    
    end
    
    # run the optimisation
    return _RBF_hypers_opt(samples,plan,kerns,rbf_opt_gens,sr,optHypRes,options) 
end


function RBFHypersResult(res,samples,kerns,variable_kernel_width::Bool,variable_dim_scaling::Bool)
    # Save the results
    nObs = length(samples)
    bestres = res.archive_output.best_candidate

    if variable_kernel_width
        x = bestres[1:nObs]     
        y = bestres[nObs+1:nObs+nObs] 
        variable_dim_scaling ? axisScale = bestres[nObs+nObs+1:end] : axisScale = 1.0

        kernind = round.(Int,_scale(y,1,length(kerns),oldMin=0,oldMax=1))
        kern = Array{ScatteredInterpolation.RadialBasisFunction,1}(undef,nObs)
        for j = 1:nObs
            kern[j] = kerns[kernind[j]](x[j])
        end
    elseif !variable_kernel_width
        x = bestres[1]
        y = bestres[2]
        variable_dim_scaling ? axisScale = bestres[3:end] : axisScale = 1.0

        kernind = round.(Int,_scale(y,1,length(kerns),oldMin=0,oldMax=1))
        kern = kerns[kernind](x)
    end


    kern_hyp_res = RBFHypersResult(
        x,                                                          #width
        kern,                                                       #kernelFunc
        axisScale,                                                  #scaling
        res.archive_output.best_fitness,                            #fitness
        )

    return kern_hyp_res
end





function _RBF_hypers_opt(   samples::Array{Float64,2},plan::Array{Float64,2},
                            kerns,rbf_opt_gens,sr,optHypRes,options)

    @unpack rippa, variable_kernel_width, variable_dim_scaling, rbf_opt_method, 
    max_rbf_width, max_scale, cond_max, rbf_dist_metric = options

    samples = vec(samples)
        
    res = bboptimize(x -> interp_obj(x,kerns,samples,plan,rippa=rippa,
            variable_kernel_width=variable_kernel_width,
            variable_dim_scaling=variable_dim_scaling,cond_max=cond_max); 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:silent, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6);
    
    bestres = res.archive_output.best_candidate;
        
    # Order the results in an Array of RBFHypersResult
    kern_hyp_res = RBFHypersResult(res,samples,kerns,
    variable_kernel_width,variable_dim_scaling)
    

    return kern_hyp_res
end




"""
    _surrogate_interpolant(optres::T,plan,samples,estimationpoints,
    oldMin,oldMax) where T <: SurrogateModelOptim.RBFHypersResult{U,Float64} where U

Evaluate an optimised interpolant at locations estimationpoints.
"""
function _surrogate_interpolant(optres::T,plan,samples,estimationpoints,
                        ) where T <: SurrogateModelOptim.RBFHypersResult{U,Float64} where U
    
    #Interpolation object based on the optimisation results
    itp = interpolate(optres.kernelFunc, plan, samples)

    #evaluate the interpolation object
    ests = similar(estimationpoints[1:1,:])
    for i = 1:size(estimationpoints,2)
        ests[i] = ScatteredInterpolation.evaluate(itp, estimationpoints[:,i])[1]
    end

    return ests
end


function _surrogate_interpolant(optres::T,samples,plan,estimationpoints
                        ) where T <: SurrogateModelOptim.RBFHypersResult{U,Array{Float64,1}} where U
    
    scaledPoints = similar(plan)
    for i = 1:size(scaledPoints,1)
        scaledPoints[i,:] = _scale(plan[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i])
    end
    
    #Interpolation object based on the optimisation results
    itp = interpolate(optres.kernelFunc, scaledPoints, samples)


    #evaluate the interpolation object
    estimationpointsScaled = similar(estimationpoints)
    for i = 1:size(estimationpointsScaled,1)
        estimationpointsScaled[i,:] = _scale(estimationpoints[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i])        
    end


    ests = similar(estimationpoints[1:1,:])
    for i = 1:size(estimationpoints,2)
        ests[i] = ScatteredInterpolation.evaluate(itp, estimationpointsScaled[:,i])[1]
    end

    return ests
end





"""
    _surrogate_interpolant(optres::T,points,observations,estimationpoints,
    oldMin,oldMax) where T <: RBFoptim_v1.HypersResult{U,Float64} where U

Evaluate an optimised interpolant at locations estimationpoints.
"""
function _surrogate_interpolant(optres::T,points,observations,estimationpoints,
    oldMin,oldMax) where T <: SurrogateModelOptim.RBFHypersResult{U,Float64} where U
    
    #Interpolation object based on the optimisation results
    itp = interpolate(optres.kernelFunc, points, observations)

    #evaluate the interpolation object
    ests = similar(estimationpoints[1:1,:])
    for i = 1:size(estimationpoints,2)
        ests[i] = ScatteredInterpolation.evaluate(itp, estimationpoints[:,i])[1]
    end

    return ests
end


function _surrogate_interpolant(optres::T,points,observations,estimationpoints,
    oldMin,oldMax) where T <: SurrogateModelOptim.RBFHypersResult{U,Array{Float64,1}} where U
    
    scaledPoints = similar(points)
    for i = 1:size(scaledPoints,1)
        scaledPoints[i,:] = _scale(points[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i],
        oldMin = oldMin[i], oldMax = oldMax[i])
    end

    
    #Interpolation object based on the optimisation results
    itp = interpolate(optres.kernelFunc, scaledPoints, observations)


    #evaluate the interpolation object
    estimationpointsScaled = similar(estimationpoints)
    for i = 1:size(estimationpointsScaled,1)
        estimationpointsScaled[i,:] = _scale(estimationpoints[i,:],-1.0*optres.scaling[i],1.0*optres.scaling[i],
        oldMin = oldMin[i], oldMax = oldMax[i])        
    end


    ests = vec(similar(estimationpoints[1:1,:]))
    for i = 1:size(estimationpoints,2)
        ests[i] = ScatteredInterpolation.evaluate(itp, estimationpointsScaled[:,i])[1]
    end

    return ests
end



#Returns anonymous function 
function surrogate_model(samples, plan, options)
    
    @unpack num_interpolants = options

    if num_interpolants == 1
        optres = _rbf_hypers_opt(samples, plan, options)

        estimationpoints->_surrogate_interpolant(
                                                optres,
                                                plan,
                                                samples',
                                                estimationpoints,
                                                minimum(plan,dims=2),
                                                maximum(plan,dims=2)
                                                )
    else
        _p_rbf_opt = function (x)
            _rbf_hypers_opt(samples, plan, options)
        end
        optres = pmap(_p_rbf_opt,1:num_interpolants')
        
        return function (estimationpoints)
            res = Array{Float64,2}(undef,num_interpolants,size(estimationpoints,2))
  
            for i = 1:length(optres)
                res[i,:]= _surrogate_interpolant(
                                                optres[i],
                                                plan,
                                                samples',
                                                estimationpoints,
                                                minimum(plan,dims=2),
                                                maximum(plan,dims=2)
                                                )
            end
            return SurrogateEstimate(res)
        end
    end
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





function MAEErrorLOO!(E, A, ests, LOOinds, interp::U,
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
        MAE = Inf
    else
        MAE = mean(abs.(E))
    end
    return MAE
end




function MAEErrorLOO!(E, A, ests, LOOinds, interp, samples, plan::T; 
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
        MAE = Inf
    else
        MAE = mean(abs.(E))
    end
    return MAE
end

