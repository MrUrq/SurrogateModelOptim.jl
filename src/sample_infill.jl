function minimum(res::SurrogateModelOptim.SurrogateEstimate{Array{Float64,2}})
    minimum(res.sm_estimate, dims = 1)
end

function maximum(res::SurrogateModelOptim.SurrogateEstimate{Array{Float64,2}})
    maximum(res.sm_estimate, dims = 1)
end

function std(res::SurrogateModelOptim.SurrogateEstimate{Array{Float64,2}})
    std(res.sm_estimate, dims = 1)
end

function median(res::SurrogateModelOptim.SurrogateEstimate{Array{Float64,2}})
    median(res.sm_estimate, dims = 1)
end

function mean(res::SurrogateModelOptim.SurrogateEstimate{Array{Float64,2}})
    mean(res.sm_estimate, dims = 1)
end

function _nearest_point(kdtree,sm_interpolant,x)
    x = permutedims(x')
    v_tmp = [x; median(sm_interpolant(x))]
    v = SVector{length(v_tmp)}(v_tmp)
    idxs, dists = knn(kdtree, v, 1, true)
    
    return dists[1]
end

function model_infill(plan,samples,sm_interpolant,options)

    @unpack rbf_opt_gens = options
    
    #KD-tree for the plan as well as the samples from the plan
    combined_plan_samples = [plan;samples]
    kdtree = KDTree(combined_plan_samples)    


    #The search takes places in the design space
    sr = vcat(extrema(plan,dims = 2)...)
    
    #The distance is calculated by picking a point in the design space, evaluating the 
    #surrogate model at that point and checking the overall distance while 
    #taking into consideration the function output. I.e constrained to the function surface
    res = bboptimize(x -> 1/_nearest_point(kdtree,sm_interpolant,x);
            Method=:de_rand_1_bin,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:silent);

    bestres = res.archive_output.best_candidate

    #Return the point with the largest found distance 
    return permutedims(bestres')
end

