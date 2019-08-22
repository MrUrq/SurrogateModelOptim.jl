"""
    _nearest_point(kdtree,sm_interpolant,x)

Helper function to calculate the distance between a sample point `x` and
all samples point in the kdtree.
"""
function _nearest_point(kdtree,sm_interpolant,x)
    x = permutedims(x') #Arrange in column vector
    v_tmp = [x; median(sm_interpolant(x))]
    v = SVector{length(v_tmp)}(v_tmp)
    _, dists = knn(kdtree, v, 1, true)
    
    return dists[1]
end

"""
    distance_infill(plan,samples,sm_interpolant)

Returns an anonymous function that calculates the distance between
sample point `x` and the closest point in `plan`. 
The function value is included in the distance calculated by using the samples
and the estimated function value.
The distance is multiplied by -1 in order to maximise the distance between sample
point `x` and all other points.
"""
function distance_infill(plan,samples,sm_interpolant)

    #KD-tree for the plan as well as the samples from the plan
    combined_plan_samples = [plan;samples]
    kdtree = KDTree(combined_plan_samples)    

    #The distance is calculated by picking a point in the design space, evaluating the 
    #surrogate model at that point and checking the overall distance while 
    #taking into consideration the function output. I.e constrained to the function surface

    return x -> -1*_nearest_point(kdtree,sm_interpolant,x)
end

"""
    minimum_infill(sm_interpolant)

Returns an anonymous function that calculates the minimum value of the surrogate
model ensemble predictions at position `x`. 
"""
function minimum_infill(sm_interpolant)
    function (x)    
        out = minimum(sm_interpolant(x))
        return out
    end
end

"""
    median_infill(sm_interpolant)

Returns an anonymous function that calculates the median value of the surrogate
model ensemble predictions at position `x`. 
"""
function median_infill(sm_interpolant)
    function (x)    
        out = median(sm_interpolant(x))
        return out
    end
end

"""
    mean_infill(sm_interpolant)

Returns an anonymous function that calculates the mean value of the surrogate
model ensemble predictions at position `x`. 
"""
function mean_infill(sm_interpolant)
    function (x)    
        out = mean(sm_interpolant(x))
        return out
    end
end

"""
    std_infill(sm_interpolant)

Returns an anonymous function that calculates the standard deviation of the surrogate
model ensemble predictions at position `x`. The std is multiplied by -1 
in order to find the sample point `x` with the largest standard deviation.
"""
function std_infill(sm_interpolant)
    function (x)
        out = -1*std(sm_interpolant(x))
        return out
    end
end

function srbf_infill(plan,samples,sm_interpolant,w)

    df = distance_infill(plan,samples,sm_interpolant)
    med = median_infill(sm_interpolant)

    function (x)
        out = med(x) + (1-w)*(1+df(x))
        return out
    end
end

function mrbf_infill(sm_interpolant,w)

    df = std_infill(sm_interpolant)
    med = median_infill(sm_interpolant)

    function (x)
        #out = med(x) + (1-w)*(1+df(x))
        out = w*med(x) + (1-w)*df(x)
        #out = med(x) + (1-w)*df(x)
        return out
    end
end

"""
    infill_objective(sm_interpolant,plan,samples,infill_funcs::Array{Symbol,1})

Returns the infill objective function based on the names supplied 
in `infill_funcs`. Choose from `[:min,:median,:mean,:dist,:std]`.
"""
function infill_objective(sm_interpolant,plan,samples,infill_funcs::Array{Symbol,1})
    #Infill function options
    min_infill_fun = minimum_infill(sm_interpolant)
    median_infill_fun = median_infill(sm_interpolant)
    mean_infill_fun = mean_infill(sm_interpolant)
    dist_infill_fun = distance_infill(plan,samples,sm_interpolant)    
    std_infill_fun = std_infill(sm_interpolant)  
    SRBF_infill_fun03 = srbf_infill(plan,samples,sm_interpolant,0.3)
    SRBF_infill_fun05 = srbf_infill(plan,samples,sm_interpolant,0.5)
    SRBF_infill_fun08 = srbf_infill(plan,samples,sm_interpolant,0.8)
    SRBF_infill_fun095 = srbf_infill(plan,samples,sm_interpolant,0.95)
    MRBF_infill_fun005 = mrbf_infill(sm_interpolant,0.05) 
    MRBF_infill_fun01 = mrbf_infill(sm_interpolant,0.1) 
    MRBF_infill_fun015 = mrbf_infill(sm_interpolant,0.15) 
    MRBF_infill_fun02 = mrbf_infill(sm_interpolant,0.2) 
    MRBF_infill_fun025 = mrbf_infill(sm_interpolant,0.25)
    MRBF_infill_fun03 = mrbf_infill(sm_interpolant,0.3)
    MRBF_infill_fun035 = mrbf_infill(sm_interpolant,0.35)
    MRBF_infill_fun04 = mrbf_infill(sm_interpolant,0.4)
    MRBF_infill_fun045 = mrbf_infill(sm_interpolant,0.45)
    MRBF_infill_fun05 = mrbf_infill(sm_interpolant,0.5)
    MRBF_infill_fun055 = mrbf_infill(sm_interpolant,0.55)
    MRBF_infill_fun06 = mrbf_infill(sm_interpolant,0.6)
    MRBF_infill_fun065 = mrbf_infill(sm_interpolant,0.65)
    MRBF_infill_fun07 = mrbf_infill(sm_interpolant,0.7)
    MRBF_infill_fun075 = mrbf_infill(sm_interpolant,0.75)
    MRBF_infill_fun08 = mrbf_infill(sm_interpolant,0.8)
    MRBF_infill_fun085 = mrbf_infill(sm_interpolant,0.85)
    MRBF_infill_fun09 = mrbf_infill(sm_interpolant,0.9)
    MRBF_infill_fun095 = mrbf_infill(sm_interpolant,0.95)


    #Get the infill objective functions
    call(f, x) = f(x)
    library = Dict(
        :min => x -> min_infill_fun(x),
        :median => x -> median_infill_fun(x),
        :mean => x -> mean_infill_fun(x),
        :dist => x -> dist_infill_fun(x),
        :std => x -> std_infill_fun(x),
        :srbf03 => x -> SRBF_infill_fun03(x),
        :srbf05 => x -> SRBF_infill_fun05(x),
        :srbf08 => x -> SRBF_infill_fun08(x),
        :srbf095 => x -> SRBF_infill_fun095(x),
        :mrbf005 => x -> MRBF_infill_fun005(x),
        :mrbf01 => x -> MRBF_infill_fun01(x),
        :mrbf015 => x -> MRBF_infill_fun015(x),
        :mrbf02 => x -> MRBF_infill_fun02(x),
        :mrbf025 => x -> MRBF_infill_fun025(x),
        :mrbf03 => x -> MRBF_infill_fun03(x),
        :mrbf035 => x -> MRBF_infill_fun035(x),
        :mrbf04 => x -> MRBF_infill_fun04(x),
        :mrbf045 => x -> MRBF_infill_fun045(x),
        :mrbf05 => x -> MRBF_infill_fun05(x),
        :mrbf055 => x -> MRBF_infill_fun055(x),
        :mrbf06 => x -> MRBF_infill_fun06(x),
        :mrbf065 => x -> MRBF_infill_fun065(x),
        :mrbf07 => x -> MRBF_infill_fun07(x),
        :mrbf075 => x -> MRBF_infill_fun075(x),
        :mrbf08 => x -> MRBF_infill_fun08(x),
        :mrbf085 => x -> MRBF_infill_fun085(x),
        :mrbf09 => x -> MRBF_infill_fun09(x),
        :mrbf095 => x -> MRBF_infill_fun095(x)
    )
    functions_to_call = Tuple([library[s] for s in infill_funcs])
    infill_obj_fun = function (x)
        call.(functions_to_call, Ref(x))
    end
    return infill_obj_fun
end

"""
    closest_index(x_val, vals) 

Based on the allowed categorical values, returns the index to the point closest
to `x_val`.
"""
function closest_index(x_val, vals)             
    ibest = first(eachindex(vals)) 
    dxbest = abs(vals[ibest]-x_val) 
    for I in eachindex(vals) 
        dx = abs(vals[I]-x_val) 
        if dx < dxbest 
            dxbest = dx 
            ibest = I 
        end 
    end 
    ibest 
end 

"""
    colinmat(mat,vec)

Check if the vector `vec` exits as a column in `mat`. Greedy implementation
potentially not fit for high speed applications of large matrices. 
"""
function colinmat(mat,vec)
    for column in eachcol(mat)
        all(column .== vec) && return true
    end
    false
end

"""
    infill_opt(search_range,infill_iterations,num_infill_points,infill_obj_fun,infill_funcs,plan,sm_interpolant,options)

Optimizes the infill points based on the supplied `infill_funcs` criteria.
Does not allow duplicated points.
"""
function infill_opt(search_range,infill_iterations,num_infill_points,infill_obj_fun,infill_funcs,plan,sm_interpolant,options)
    # Try generating pareto optimal infill points. Wrapped in try block due to 
    # the algorithm sometimes failing to generate a result.
    infill_incomplete = true
    j = 0
    res_bboptim = nothing
    while infill_incomplete && j <= 100
        try 
            res_bboptim = bboptimize(infill_obj_fun; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{length(infill_funcs[1:num_infill_points])}(is_minimizing=true),
                    SearchRange=search_range, ϵ=0.000001,
                    MaxFuncEvals=infill_iterations,
                    MaxStepsWithoutProgress=infill_iterations,TraceMode=:silent); 
            infill_incomplete = false
        catch ex
            j += 1 
            (j >= 100) && rethrow(ex)  #Show error if consistently failing
        end
    end
    
    # Pick the infill points 
    infill_plan = Array{Float64,2}(undef,size(plan,1),0)
    infill_type = Array{Symbol,1}()
    infill_prediction = Array{Float64,1}()

    # Take the num_infill_points number of infill points from the cycled list.
    infill_obj_funs = Iterators.take(Base.Iterators.cycle(1:length(infill_funcs[1:num_infill_points])), num_infill_points)
    
    for i in infill_obj_funs
        pf = pareto_frontier(res_bboptim)
        best_obj1, idx_obj1 = findmin(map(elm -> fitness(elm)[i], pf))
        bo1_solution = BlackBoxOptim.params(pf[idx_obj1]) # get the solution candidate itself... 

        # Add the infill point if it does not exist in the plan or infill_plan
        v = copy(permutedims(bo1_solution'))
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,infill_funcs[i])
        end
    end

    # Cycle through the objective functions and update options with new list.
    ifuncs=circshift(infill_funcs,-num_infill_points)
    circshift_perm=[findfirst(isequal(x),ifuncs) for x in infill_funcs]
    options = Options(options; infill_funcs=ifuncs) 

    return infill_plan,infill_type,infill_prediction,res_bboptim,options
end

"""
    infill_add(sm_interpolant,samples,plan,infill_prediction,search_range,infill_plan,infill_type,num_infill_points,infill_obj_fun,infill_funcs,res_bboptim)

Adds additional sample points if the `infill_opt` can't supply the number of 
requested infill points. First a pareto optimal point based on the infill 
functions is selected. The pareto point furthest away from all existing points 
is selected. If this fails the last points are added as randomly selected points
to ensure that the requested number of infill points is added. Does not allow
duplicated points.
"""
function infill_add(sm_interpolant,samples,plan,infill_prediction,search_range,infill_plan,infill_type,num_infill_points,infill_obj_fun,infill_funcs,res_bboptim)
    
    for i = 1:(num_infill_points-size(infill_plan,2)) 
        infill_samples = Array{Float64,2}(undef,1,size(infill_plan,2))
        for i = 1:size(infill_plan,2)
            infill_samples[i] = median(sm_interpolant(infill_plan[:,i]))
        end

        # Find the point in the pareto front which is located furthest away from previous samples
        dist_infill_fun = distance_infill([plan infill_plan],[samples infill_samples],sm_interpolant)  
        dist_obj = Inf
        par_f = pareto_frontier(res_bboptim)
        best_add_infill_solution = BlackBoxOptim.params(par_f[1])
        for pf in par_f
            bo1_solution = BlackBoxOptim.params(pf)
            if (cur_val = dist_infill_fun(bo1_solution)) < dist_obj
                dist_obj = cur_val
                best_add_infill_solution = copy(bo1_solution)                
            end
        end
        
        # Add the infill point if it does not exist in the plan or infill_plan
        v = copy(permutedims(best_add_infill_solution'))
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,:pareto)
        end       
    end
    
    # Add random points if there aren't enough in the pareto front
    while (num_infill_points-size(infill_plan,2)) != 0
        
        v = Array{Float64,2}(undef,size(infill_plan,1),1)
        for i = 1:length(v)
            v[i] = scale(rand(),search_range[i][1],search_range[i][2],old_min=0,old_max=1)
        end
        # Add the infill point if it does not exist in the plan or infill_plan
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,:rand)
        end
    end
    
    return infill_plan,infill_type,infill_prediction
end

print_infill_head() = println("Finding new infill samples")

function print_infill_tail(infill_type,num_infill_points,infill_plan,infill_prediction)
    print("    ")
    for i = (length(infill_type)-num_infill_points+1):length(infill_type)
        print(@sprintf("%-15s",infill_type[i]))
    end
    print("\n")
    println("    ---------------------------------------------------------------")
    for j = 1:size(infill_plan,1)
        print("    ")
        for i = 1:size(infill_plan,2)
            print(@sprintf("%-15.7g",infill_plan[j,i]))
        end
        print("\n")
    end
    println("    ---------------------------------------------------------------")
    _, min_loc = findmin(infill_prediction)
    print("    ")
    for i = 1:size(infill_plan,2)
        if i == min_loc
            printstyled(@sprintf("%-15.7g",infill_prediction[i]); bold=true)
        else
            printstyled(@sprintf("%-15.7g",infill_prediction[i]))
        end
    end
    println("\t prediction\n")
end