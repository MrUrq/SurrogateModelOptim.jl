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
        end
    end
    return Tuple(output)
end

function extract_bboptim_hypers(bboptim_fcall_vector,plan,kerns,variable_kernel_width,variable_dim_scaling,smooth)

    n_dims, n_samples = size(plan)
    
    # Kernel width and type length
    variable_kernel_width ? n_kerns = n_samples : n_kerns = 1    
    
    # Dimensional scaling length
    variable_dim_scaling ? n_dim_scales = n_dims : n_dim_scales = 0

    # Ridge regression smoothing length
    (smooth == :variable) && (n_smooth = n_samples)
    (smooth == :single)   && (n_smooth = 1)
    (smooth == false)     && (n_smooth = 0)

    width, kernel_float, scaling, smooth = extract_vector_range(n_kerns,
                                                                n_kerns,
                                                                n_dim_scales,
                                                                n_smooth)
    
    kern_ind = round.(Int,_scale(kernel_float,1,length(kerns),old_min=0,old_max=1))

    return width, kern_ind, scaling, smooth
end
