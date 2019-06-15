"""
    _nearest_point(kdtree,sm_interpolant,x)

Helper function to calculate the distance between a sample point `x` and
all samples point in the kdtree.
"""
function _nearest_point(kdtree,sm_interpolant,x)
    x = permutedims(x')
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

    #Get the infill objective functions
    call(f, x) = f(x)
    library = Dict(
        :min => x -> min_infill_fun(x),
        :median => x -> median_infill_fun(x),
        :mean => x -> mean_infill_fun(x),
        :dist => x -> dist_infill_fun(x),
        :std => x -> std_infill_fun(x)
    )
    functions_to_call = Tuple([library[s] for s in infill_funcs])
    infill_obj_fun = function (x)
        call.(functions_to_call, Ref(x))
    end
    return infill_obj_fun
end