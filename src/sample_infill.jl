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

function _nearest_point(kdtree,res,sm_interpolant,x)
    x = permutedims(x')
    v_tmp = [x; median(sm_interpolant(res,x))]
    v = SVector{length(v_tmp)}(v_tmp)
    idxs, dists = knn(kdtree, v, 1, true)
    
    return minimum((1e10,dists[1]))    
end

function distance_infill(plan,samples,res,sm_interpolant)

    #KD-tree for the plan as well as the samples from the plan
    combined_plan_samples = [plan;samples]
    kdtree = KDTree(combined_plan_samples)    

    #The distance is calculated by picking a point in the design space, evaluating the 
    #surrogate model at that point and checking the overall distance while 
    #taking into consideration the function output. I.e constrained to the function surface

    return x -> 1/_nearest_point(kdtree,res,sm_interpolant,x)
end

function minimum_infill(res,sm_interpolant)
    function (x)    
        out = median(sm_interpolant(res,x))[1]
        return minimum((1e10,out))
    end
end

function std_infill(res,sm_interpolant)
    function (x)
        out = 1/(std(sm_interpolant(res,x)))[1]
        return minimum((1e10,out))
    end
end

function min_std_infill(c,res,sm_interpolant)
    function (x)
        ret = sm_interpolant(res,x)
        out = median(ret)[1]-c*std(ret)[1]

        return minimum((1e10,out))
    end
end

function min_std_zscore_infill(c,res,sm_interpolant)
    function (x)
        y = sm_interpolant(res,x).sm_estimate

        y_filt = y[findall((x -> (x < 1) & (x > -1)), zscore(y))]
        out = median(y)[1]-c*std(y_filt)[1]

        return minimum((1e10,out))
    end
end


function model_infill(plan,samples,sm_interpolant,criteria,options)

    @unpack rbf_opt_gens, num_interpolants, num_infill_points,
            trace, infill_funcs = options

    if trace
        println("Finding new infill samples ...")
    end

    res_alloc = Array{Float64,2}(undef,num_interpolants,1)
    
    #Infill function options
    min_infill_fun = minimum_infill(res_alloc,sm_interpolant)
    min_std_infill_fun = min_std_infill(2,res_alloc,sm_interpolant)
    dist_infill_fun = distance_infill(plan,samples,res_alloc,sm_interpolant)    
    std_infill_fun = std_infill(res_alloc,sm_interpolant)     
    min_std_infill_zscore_fun = min_std_zscore_infill(1,res_alloc,sm_interpolant) 
    
    #The search takes places in the design space
    sr = vcat(extrema(plan,dims = 2)...)

    #Get the infill objective functions
    call(f, x) = f(x)
    library = Dict(
        :min => x -> min_infill_fun(x),
        :min_2std => x -> min_std_infill_fun(x),
        :dist => x -> dist_infill_fun(x),
        :std => x -> std_infill_fun(x),
        :min_std_z => x -> min_std_infill_zscore_fun(x),
    )
    functions_to_call = Tuple([library[s] for s in infill_funcs])
    infill_obj_fun = function (x)
        call.(functions_to_call, Ref(x))
    end
    
    num_infill_obj_funs = length(infill_obj_fun(plan[:,1]))

    # Try generating pareto optimal infill points. Wrapped in try block due to 
    # the algorithm sometimes failing to generate a result.
    infill_incomplete = true
    j = 0
    res_bboptim = nothing
    while infill_incomplete && j < 51
        try 
            res_bboptim = bboptimize(infill_obj_fun; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{num_infill_obj_funs}(is_minimizing=true),
                    SearchRange=sr, ϵ=0.00001,
                    MaxFuncEvals=50000, TraceMode=:silent);
            infill_incomplete = false
        catch
            j += 1 
        end
    end
    

    # Pick the infill points 
    infill_plan = Array{Float64,2}(undef,size(plan,1),0)
    infill_type = Array{Symbol,1}()
    infill_prediction = Array{Float64,1}()

    # Cycle through the objective functions if the number of infill points is 
    # less than the number of objective functions.
    if num_infill_points < num_infill_obj_funs
        
        infill_obj_funs = Array{Int64}(undef,num_infill_points)
        for i = 1:num_infill_points
            infill_obj_funs[i] = first(Iterators.drop(Iterators.cycle(1:num_infill_obj_funs),criteria+num_infill_obj_funs-2+i))
        end

        # Find the next cyclic infill objective function to use
        criteria = infill_obj_funs[end]+1
    else
        infill_obj_funs = 1:num_infill_obj_funs
    end
    
    for i in infill_obj_funs
        pf = pareto_frontier(res_bboptim)
        best_obj1, idx_obj1 = findmin(map(elm -> fitness(elm)[i], pf))
        bo1_solution = params(pf[idx_obj1]) # get the solution candidate itself... 

        # Do not add duplicate points
        v = copy(permutedims(bo1_solution'))        
        if !any(c->view(infill_plan,:,c)==v,1:size(infill_plan,2))
            push!(infill_prediction,median(sm_interpolant(res_alloc,bo1_solution))[1])
            infill_plan = [infill_plan v]            
            push!(infill_type,infill_funcs[i])
        end
    end

    # Add additional points to fill up the desired number of infill points   
    while (num_infill_points-size(infill_plan,2)) != 0
        infill_samples = Array{Float64,2}(undef,1,size(infill_plan,2))
        for i = 1:size(infill_plan,2)
            infill_samples[i] = median(sm_interpolant(res_alloc,infill_plan[:,i]))[1]
        end

        # Find the point in the pareto front which is located furthest away from previous samples
        dist_infill_fun = distance_infill([plan infill_plan],[samples infill_samples],res_alloc,sm_interpolant)  
        dist_obj = Inf
        par_f = pareto_frontier(res_bboptim)
        best_add_infill_solution = params(par_f[1])
        for pf in par_f
            bo1_solution = params(pf) 
            if (cur_val = dist_infill_fun(bo1_solution)) < dist_obj
                dist_obj = cur_val
                best_add_infill_solution = copy(bo1_solution)                
            end
        end
        
        v = copy(permutedims(best_add_infill_solution'))    
        if !any(c->view(infill_plan,:,c)==v,1:size(infill_plan,2))
            push!(infill_prediction,median(sm_interpolant(res_alloc,best_add_infill_solution))[1])
            infill_plan = [infill_plan v]
            push!(infill_type,:pareto)
        end       
    end

    if trace
        println("Infill samples:")
        println("----------------------------------------------")
        for i = (length(infill_type)-num_infill_points+1):length(infill_type)
            print(@sprintf("%-10s",infill_type[i]))
        end
        print("\n")
        println("----------------------------------------------")
        for j = 1:size(infill_plan,1)
            for i = 1:size(infill_plan,2)
                print(@sprintf("%-10.04f",infill_plan[j,i]))
            end
            print("\n")
        end
        println("----------------------------------------------")
        for i = 1:size(infill_plan,2)
            printstyled(@sprintf("%-10.05f",infill_prediction[i]); color=:light_green)
        end
        print("\t prediction\n")
        println("----------------------------------------------")
    end

    return infill_plan, criteria, infill_type, infill_prediction
end

