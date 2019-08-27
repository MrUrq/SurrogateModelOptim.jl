"""
    function _rippa(A,a)

Estimate the Leave-One-Out (LOO) errors using rippas method. Complexity of 𝛰(3)
compared to calculating the exact LOO at a cost of 𝛰(4). `A` is the RBF matrix
and `a` is the weights of the RBF.
"""
function _rippa(A, a)
    N = size(A, 1)           #Number of Leave-One-Out (LOO) errors
    xₖ = similar(A)
    E = Array{Float64}(undef,N)

    _rippa!(E,xₖ,A,a)
end


"""
    function _rippa!(E,xₖ,A,a)

Same as _rippa but inplace.
"""
function _rippa!(E,xₖ,A,a)
    N = size(A, 1)           #Number of Leave-One-Out (LOO) errors
    xₖ = A\I
    @inbounds for k = 1:N
        E[k] = a[k] / xₖ[k,k]   #Estimated error for k-th subset
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
    cond_max::Float64=1e6, cond_check::Bool=false, rippa::Bool=false, rbf_dist_metric = Euclidean())

    #initiate arrays
    N = length(samples)
    E = Array{Float64}(undef,N)       #Error
    A = zeros(Float64, N, N)    #RBF matrix A
    ests = Array{Float64}(undef,N)    #Function estimate based on LOO
    LOOinds = Array{Int}(undef,N - 1, N)
    for i = 1:N
        LOOinds[:,i] = filter(x -> x != i, 1:N) # Get the leave one out sub indices
    end
    xₖ = similar(A)

    RMSErrorLOO!(E, A, ests, xₖ, LOOinds, interp, samples, plan, smooth;
    cond_max=cond_max, rippa=rippa, rbf_dist_metric = rbf_dist_metric)
end


function RMSErrorLOO!(E, A, ests, xₖ, LOOinds, interp::U,
 samples, plan::T, smooth; cond_max::Float64=1e6, cond_check::Bool=false,
  rippa::Bool=false, rbf_dist_metric = Euclidean()) where T <: AbstractArray where U <: AbstractArray

    N = length(samples)
    RMSE = Inf

    #perform Leave-One-Out estimation and calculate error
    if rippa
        @assert typeof(interp) <: Vector{ScatteredInterpolation.RadialBasisFunction} "Rippas
              algorithm only available for Radial Basis Functions"

        #interpolation object trained on the entire dataset
        itp, A = interpolate(interp, plan,
            samples, returnRBFmatrix=true, metric = rbf_dist_metric, smooth=smooth)

        #check the conditioning of matrix A for RBFs
        if cond_check && (cond(A) > cond_max)
            E .= Inf
        else
            #RBF error estimation based on Rippas algorithm
            E = _rippa!(E, xₖ, A, itp.w)
            RMSE = sqrt(mean(E.^2))
        end
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

            #check the conditioning of matrix A for RBFs
            if cond_check && (cond(A) > cond_max)
                E .= Inf
            else
                #evaluate the interpolation object in the LOO position
                ests[i] = ScatteredInterpolation.evaluate(itp, plan[:,i])[1]

                E[i] = samples[i] - ests[i]
                RMSE = sqrt(mean(E.^2))
            end
        end
    end

    return RMSE
end

function RMSErrorLOO!(E, A, ests, xₖ, LOOinds, interp, samples, plan::T, smooth;
 cond_max::Float64=1e6, cond_check::Bool=false, rippa::Bool=false, rbf_dist_metric = Euclidean()) where T <: AbstractArray

    N = length(samples)
    RMSE = Inf
    
    #perform Leave-One-Out estimation and calculate error
    if rippa
        @assert typeof(interp) <: ScatteredInterpolation.RadialBasisFunction "Rippas algorithm only available for Radial Basis Functions"

        #interpolation object trained on the entire dataset
        itp, A = interpolate(interp, plan,
         samples, returnRBFmatrix=true, metric = rbf_dist_metric, smooth=smooth)

        #check the conditioning of matrix A for RBFs
        if cond_check && (cond(A) > cond_max)
            E .= Inf
        else
            #RBF error estimation based on Rippas algorithm
            E = _rippa!(E, xₖ, A, itp.w)
            RMSE = sqrt(mean(E.^2))
        end
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
            
            #check the conditioning of matrix A for RBFs
            if cond_check && (cond(A) > cond_max)
                E .= Inf
            else
               #evaluate the interpolation object in the LOO position
                ests[i] = ScatteredInterpolation.evaluate(itp, plan[:,i])[1]

                E[i] = samples[i] - ests[i]
                RMSE = sqrt(mean(E.^2))
            end
        end
    end

    return RMSE    
end

"""
    function scale(old_x::T,new_min,new_max;old_min=minimum(old_x),
    old_max=maximum(old_x)) where T <: Real

Scale a scalar to match a new range.
"""
function scale(old_x::T,new_min,new_max;old_min=minimum(old_x),old_max=maximum(old_x)) where T <: Real

    old_x = (((old_x - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
end

"""
    function scale(old_x::Array{T,1},new_min,new_max;old_min=minimum(old_x),
    old_max=maximum(old_x)) where T <: Real

Scale a vector to match a new range.
"""
function scale(old_x::Array{T,1},new_min,new_max;old_min=minimum(old_x),old_max=maximum(old_x)) where T <: Real

    newX = (((old_x .- old_min) .* (new_max .- new_min)) ./ (old_max .- old_min)) .+ new_min
end

"""
    function scale(old_x::Array{T,2},direction::Int,new_min,new_max;
    old_min=minimum(old_x,direction)::Array{T,2},old_max=maximum(old_x,direction)::Array{T,2}) where T <: Real 

Scale a 2D matrix to match a new range along the specified `direction`.
"""
function scale(old_x::Array{T,2},direction::Int,new_min,new_max;old_min=minimum(old_x,dims=direction)::Array{T,2},old_max=maximum(old_x,dims=direction)::Array{T,2}) where T <: Real 

    newX = mapslices(x -> scale(x,new_min,new_max;old_min=old_min,old_max=old_max), old_x, dims=direction)
end

scale(x::Missing,min_val,max_val;old_min,old_max) = missing

"""
    preprocess_point(points,scaling;base_scale::Array{Float64,2})

Scales the input points in order to evaluate the adaptively scaled RBF correctly.
"""
function preprocess_point(points,scaling;base_scale::Array{Float64,2})

    min_max = extrema(base_scale, dims=2)

    preprocessed_point = similar(points)
    for i = 1:size(preprocessed_point,1)
        preprocessed_point[i,:] = scale(points[i,:],-1.0*scaling[i],1.0*scaling[i],
        old_min = min_max[i][1], old_max = min_max[i][2])
    end
    return preprocessed_point
end

"""
    preprocess_point!(preprocessed_point,points,scaling;old_min::Array{Float64,2},old_max::Array{Float64,2})

Scales the input points in order to evaluate the adaptively scaled RBF correctly.
"""
function preprocess_point!(preprocessed_point,points,scaling;old_min::Array{Float64,2},old_max::Array{Float64,2})

    for i = 1:size(preprocessed_point,1)
        preprocessed_point[i,:] = scale(points[i,:],-1.0*scaling[i],1.0*scaling[i],
        old_min = old_min[i], old_max = old_max[i])
    end
    return preprocessed_point
end

"""
    preprocess_point!(preprocessed_point,points,scaling::Bool;old_min::Array{Float64,2},old_max::Array{Float64,2}) where U

Scales the input points in order to evaluate the adaptively scaled RBF correctly.
"""
function preprocess_point!(preprocessed_point,points,scaling::Bool;old_min::Array{Float64,2},old_max::Array{Float64,2}) where U
    preprocessed_point = points
end

"""
    preprocess_point(points,scaling::Bool;base_scale::Array{Float64,2}) where U

Scales the input points in order to evaluate the adaptively scaled RBF correctly.
"""
function preprocess_point(points,scaling::Bool;base_scale::Array{Float64,2}) where U
    points    
end

"""
    surrogate_evaluate(preprocessed_est_point,estimation_point,itp,optres,old_min,old_max)

Processes the input to scale correctly based on the options used.
"""
function surrogate_evaluate(preprocessed_est_point,estimation_point,itp,optres,old_min,old_max)
    preprocessed_est_point = preprocess_point!( preprocessed_est_point,
                                                estimation_point, optres.scaling;
                                                old_min=old_min,
                                                old_max=old_max,
                                                )
    
    # Evaluate the interpolation object
    res = ScatteredInterpolation.evaluate(itp, preprocessed_est_point)[1]
end

"""
    construct_search_range(plan::Array{Float64,2}, variable_kernel_width,
    min_rbf_width, max_rbf_width, variable_dim_scaling,
    min_scale, max_scale, smooth, max_smooth)
    
Creates array of tuples for the optimization of the RBF hyperparameters.
"""
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

"""
    create_sr(vargs::NamedTuple{(:min_range, :max_range, :n_times),Tuple{Float64,Float64,Int64}}...)

Facilitate the creation of array of tuples for the optimization of the RBF hyperparameters.
"""
function create_sr(vargs::NamedTuple{(:min_range, :max_range, :n_times),Tuple{Float64,Float64,Int64}}...)
    sr = Array{Tuple{Float64,Float64},1}()
    for varg in vargs
        for i = 1:varg[3]
            push!(sr,(varg[1],varg[2]))
        end
    end        
    return sr
end

"""
    extract_vector_range(vargs::Int64...)   

Supply length of subsets contained in a vector and receive a tuple containing
all the ranges needed to extract each subset in the vector.
"""
function extract_vector_range(vargs::Int...)
    output = Array{Union{UnitRange{Int}, Int, Bool},1}()

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

"""
    extract_bboptim_hypers(bboptim_fcall_vector,plan,kerns,
    variable_kernel_width,variable_dim_scaling,
    smooth,smooth_user)   

Get the optimized RBF hyperparameters in useful format for further use.
"""
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
    kern_ind = round.(Int,scale(kernel_float,1,length(kerns),old_min=0,old_max=1))
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

"""
    interp_obj(inpt::Vector{Float64}, kerns, samples,
    plan::Array{Float64,2}; rippa::Bool = false,
    variable_kernel_width::Bool = true, variable_dim_scaling::Bool = true,
    smooth = false, cond_max=cond_max, rbf_dist_metric = Distances.Euclidean(),
    smooth_user::Float64 = 0.0)

Objective function for optimization of interpolation kernel function and width.
"""
function interp_obj(inpt::Vector{Float64}, E, A, ests, e, LOOinds, kerns, samples,
        plan::Array{Float64,2}; rippa::Bool = false,
        variable_kernel_width::Bool = true, variable_dim_scaling::Bool = true,
        smooth::Union{Bool, Symbol} = false, cond_max::Float64=Inf,  cond_check::Bool=false, rbf_dist_metric = Distances.Euclidean(),
        smooth_user::Float64 = 0.0)

    
    kern, scaling, smooth = extract_bboptim_hypers( inpt,plan,kerns,variable_kernel_width,
                                                    variable_dim_scaling,smooth,
                                                    smooth_user)

    #Preprocess the plan based on the settings used.
    preprocessed_plan = preprocess_point(plan,scaling,base_scale=plan)
    
    #Wrapped in try catch-block to catch failure to solve linear eq. system
    E = try
        E = RMSErrorLOO!(E, A, ests, e, LOOinds, kern, samples, preprocessed_plan, smooth; rippa = rippa,
        cond_max = cond_max, rbf_dist_metric = rbf_dist_metric, cond_check=cond_check)
    catch 
        E = Inf
    end

    return E
end

"""
    rbf_hypers_opt(samples_org::Array{Float64,2}, plan::Array{Float64,2}, options::Options)

Optimization function of Radial Basis Function kernel and width.
"""
function rbf_hypers_opt(samples_org::Array{Float64,2}, plan::Array{Float64,2}, options::Options)
    
    @unpack rippa, variable_kernel_width, variable_dim_scaling, rbf_opt_method, 
            min_rbf_width, max_rbf_width, min_scale, max_scale, cond_max,
            rbf_dist_metric, rbf_opt_gens, kerns, rbf_opt_pop,
            smooth, max_smooth, smooth_user, cond_check = options

    # Sample order for ScatteredInterpolation
    samples = vec(samples_org)














    # Create the hyperparameter search range based on the input options 
    sr = construct_search_range(plan, false,
                                min_rbf_width, max_rbf_width, true,
                                min_scale, max_scale, smooth, max_smooth)

                        
    # Pre-allocate for inplace functions    
    N = length(samples) #Number of Leave-One-Out (LOO) errors
    e = [[(i == j) ? true : false for i in 1:N] for j in 1:N]
    E = Array{Float64}(undef,N)
    A = zeros(Float64, N, N)    #RBF matrix A
    ests = Array{Float64}(undef,N)    #Function estimate based on LOO
    LOOinds = Array{Int}(undef,N - 1, N)
    for i = 1:N
        LOOinds[:,i] = filter(x -> x != i, 1:N) # Get the leave one out sub indices
    end
    
    # RBF hyperparameter objective function
    itp_obj = function (x)
        interp_obj(x,E,A,ests,e,LOOinds,kerns,samples,plan; 
                rippa=rippa, variable_kernel_width=false,
                variable_dim_scaling=variable_dim_scaling, smooth=smooth,
                cond_max=cond_max,rbf_dist_metric=rbf_dist_metric,cond_check=cond_check)
    end

    # Optimize the interpolant hyperparameters
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:silent, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6,
            PopulationSize = rbf_opt_pop,
            MaxStepsWithoutProgress=rbf_opt_gens,
            MaxNumStepsWithoutFuncEvals=rbf_opt_gens,
            );
        
    kern, scaling, smoothing = extract_bboptim_hypers( res.archive_output.best_candidate,
                                                    plan,kerns,false,
                                                    variable_dim_scaling,smooth,smooth_user)
    #@show best_fitness(res)




















    # Create the hyperparameter search range based on the input options 
    sr = construct_search_range(plan, variable_kernel_width,
                                min_rbf_width, max_rbf_width, variable_dim_scaling,
                                min_scale, max_scale, smooth, max_smooth)

            
    # Pre-allocate for inplace functions    
    N = length(samples) #Number of Leave-One-Out (LOO) errors
    e = [[(i == j) ? true : false for i in 1:N] for j in 1:N]
    E = Array{Float64}(undef,N)
    A = zeros(Float64, N, N)    #RBF matrix A
    ests = Array{Float64}(undef,N)    #Function estimate based on LOO
    LOOinds = Array{Int}(undef,N - 1, N)
    for i = 1:N
        LOOinds[:,i] = filter(x -> x != i, 1:N) # Get the leave one out sub indices
    end
    
    # RBF hyperparameter objective function
    itp_obj = function (x)
        interp_obj(x,E,A,ests,e,LOOinds,kerns,samples,plan; 
                rippa=rippa, variable_kernel_width=variable_kernel_width,
                variable_dim_scaling=variable_dim_scaling, smooth=smooth,
                cond_max=cond_max,rbf_dist_metric=rbf_dist_metric,cond_check=cond_check)
    end


    x0 = res.archive_output.best_candidate
    #kern    scaling    smooth
    x1=[[x0[1] for _ in 1:size(plan,2)]; [x0[2] for _ in 1:size(plan,2)]; x0[3:end]]
    nvals=length(sr)
    @assert length(x1)==length(sr)
    pop = [x1[i] for i in 1:nvals, _ in 1:rbf_opt_pop÷2]

    pop_rand = similar(pop)
    for i in 1:size(pop_rand,2)
        pop_rand[:,i] = BlackBoxOptim.Utils.latin_hypercube_sampling(first.(sr),last.(sr),1)
    end
    pop = [pop pop_rand]


    # x0 = res.archive_output.best_candidate
    # #kern    scaling    smooth
    # x1=[[x0[1] for _ in 1:size(plan,2)]; [x0[2] for _ in 1:size(plan,2)]; x0[3:end]]
    # nvals=length(sr)
    # @assert length(x1)==length(sr)
    # pop = [x1[i] for i in 1:nvals, _ in 1:1]

    # pop_rand = Array{Float64,2}(undef,length(pop),rbf_opt_pop-1)
    # for i in 1:size(pop_rand,2)
    #     pop_rand[:,i] = BlackBoxOptim.Utils.latin_hypercube_sampling(first.(sr),last.(sr),1)
    # end
    # pop = [pop pop_rand]


    # Optimize the interpolant hyperparameters
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:silent, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6,
            PopulationSize = rbf_opt_pop,
            MaxStepsWithoutProgress=rbf_opt_gens,
            MaxNumStepsWithoutFuncEvals=rbf_opt_gens,
            Population=pop
            );
    #@show best_fitness(res)
    kern, scaling, smoothing = extract_bboptim_hypers( res.archive_output.best_candidate,
                                                    plan,kerns,variable_kernel_width,
                                                    variable_dim_scaling,smooth,smooth_user)
    
    # Return the optimized hyperparameters in the correct type
    return RBFHypers(kern, scaling, smoothing), best_fitness(res)
end