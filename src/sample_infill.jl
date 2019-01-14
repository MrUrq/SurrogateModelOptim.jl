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

function distance_infill(plan,samples,sm_interpolant,options)

    #KD-tree for the plan as well as the samples from the plan
    combined_plan_samples = [plan;samples]
    kdtree = KDTree(combined_plan_samples)    

    #The distance is calculated by picking a point in the design space, evaluating the 
    #surrogate model at that point and checking the overall distance while 
    #taking into consideration the function output. I.e constrained to the function surface

    return x -> 1/_nearest_point(kdtree,sm_interpolant,x)
end

function minimum_infill(sm_interpolant)
    x->median(sm_interpolant(x))[1]
end

function std_infill(sm_interpolant)
    x->1/(std(sm_interpolant(x)))[1]
end

function min_std_infill(c,sm_interpolant)
    x->median(sm_interpolant(x))[1]-c*std(sm_interpolant(x))[1]
end

function min_std_zscore_infill(c,sm_interpolant)
    function (x)
        y = sm_interpolant(x).sm_estimate

        y_filt = y[findall((x -> (x < 1) & (x > -1)), zscore(y))]
        median(y)[1]-c*std(y_filt)[1]
    end
end


function model_infill(plan,samples,sm_interpolant,options)

    @unpack rbf_opt_gens = options


    dist_infill_fun = distance_infill(plan,samples,sm_interpolant,options)

    min_infill_fun = minimum_infill(sm_interpolant)

    std_infill_fun = std_infill(sm_interpolant) 

    min_std_infill_fun = min_std_infill(2,sm_interpolant)

    min_std_infill_zscore_fun = min_std_zscore_infill(1,sm_interpolant)

    #The search takes places in the design space
    sr = vcat(extrema(plan,dims = 2)...)
    

               
    
    
    
    function fitness_all(x)
        return (minimum((1e10,dist_infill_fun(x))),
                minimum((1e10,min_infill_fun(x))),
                minimum((1e10,std_infill_fun(x))),
                minimum((1e10,min_std_infill_fun(x))))
    end
    
    infill_incomplete = true
    while infill_incomplete
        try 
            global res = bboptimize(fitness_all; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{4}(is_minimizing=true),
                    SearchRange=sr, ϵ=0.001,
                    MaxFuncEvals=100000, TraceInterval=1.0, TraceMode=:silent);
            infill_incomplete = false
        catch
        end
    end
    
    for i = 1:4
        pf = pareto_frontier(res)
        best_obj1, idx_obj1 = findmin(map(elm -> fitness(elm)[i], pf))
        bo1_solution = params(pf[idx_obj1]) # get the solution candidate itself... 
        @show bo1_solution
        @show best_obj1
        plan = [plan deepcopy(permutedims(bo1_solution'))]
    end
    return plan   




    # fun = [x->minimum((1e10,min_infill_fun(x))),
    #        x->minimum((1e10,min_std_infill_zscore_fun(x))),
    #        x->minimum((1e10,dist_infill_fun(x)))]

    # for i = 1:length(fun)
    #     res = bboptimize(fun[i];
    #             Method=:de_rand_1_bin,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
    #             TraceMode=:silent);
    #     @show bestres = res.archive_output.best_candidate
    #     @show bestfit = res.archive_output.best_fitness
        
    #     #Return the point with the best function value
    #     plan = [plan permutedims(bestres')]        
    # end
    return plan   
end

