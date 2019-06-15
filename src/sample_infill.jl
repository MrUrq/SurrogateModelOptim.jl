"""
    _nearest_point(kdtree,sm_interpolant,x)

Compute the Bar index between `x` and `y`. If `y` is missing, compute
the Bar index between all pairs of columns of `x`.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function model_infill(sr::Array{Tuple{Float64,Float64},1},plan::AbstractArray{T,2},
        samples::AbstractArray{T,2},sm_interpolant,criteria,options) where T

    @unpack rbf_opt_gens, num_interpolants, num_infill_points,
            trace, infill_funcs, infill_iterations = options

    if trace
        println("Finding new infill samples ...")
    end
    
    infill_obj_fun = infill_objective(sm_interpolant,plan,samples,infill_funcs)
    num_infill_obj_funs = length(infill_obj_fun(plan[:,1]))

    # Try generating pareto optimal infill points. Wrapped in try block due to 
    # the algorithm sometimes failing to generate a result.
    infill_incomplete = true
    j = 0
    res_bboptim = nothing
    
    while infill_incomplete && j <= 50
        try 
            res_bboptim = bboptimize(infill_obj_fun; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{num_infill_obj_funs}(is_minimizing=true),
                    SearchRange=sr, ϵ=0.00001,
                    MaxFuncEvals=infill_iterations,
                    MaxStepsWithoutProgress=20_000,TraceMode=:silent); 
            infill_incomplete = false
        catch ex
            j += 1 
            (j >= 50) && rethrow(ex)  #Show error if consistently failing
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
        bo1_solution = BlackBoxOptim.params(pf[idx_obj1]) # get the solution candidate itself... 

        # Add the infill point if it does not exist in the plan or infill_plan
        v = copy(permutedims(bo1_solution'))
        # Add the infill point if it does not exist in the plan or infill_plan
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,infill_funcs[i])
        end
    end

    # Add additional points to fill up the desired number of infill points   
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
        # Add the infill point if it does not exist in the plan or infill_plan
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,:pareto)
        end       
    end
    
    # Add random points if there aren't enough in the pareto front
    while (num_infill_points-size(infill_plan,2)) != 0

        # Add the infill point if it does not exist in the plan or infill_plan
        v = Array{Float64,2}(undef,size(infill_plan,1),1)
        for i = 1:length(v)
            v[i] = _scale(rand(),sr[i][1],sr[i][2],old_min=0,old_max=1)
        end

        # Add the infill point if it does not exist in the plan or infill_plan
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,:rand)
        end
    end

    if trace
        println("Infill samples:")
        println("---------------------------------------------------------------")
        for i = (length(infill_type)-num_infill_points+1):length(infill_type)
            print(@sprintf("%-15s",infill_type[i]))
        end
        print("\n")
        println("---------------------------------------------------------------")
        for j = 1:size(infill_plan,1)
            for i = 1:size(infill_plan,2)
                print(@sprintf("%-15.7g",infill_plan[j,i]))
            end
            print("\n")
        end
        println("---------------------------------------------------------------")
        _, min_loc = findmin(infill_prediction)
        for i = 1:size(infill_plan,2)
            if i == min_loc
                printstyled(@sprintf("%-15.7g",infill_prediction[i]); color=:light_green, bold=true)
            else
                printstyled(@sprintf("%-15.7g",infill_prediction[i]))
            end
        end
        print("\t prediction\n")
        println("---------------------------------------------------------------")
    end

    return infill_plan, criteria, infill_type, infill_prediction
end

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

function colinmat(mat,vec)
    for column in eachcol(mat)
        all(column .== vec) && return true
    end
    false
end








































function model_infill(sr,plan,samples,sm_interpolant,vals,criteria,options)

    @unpack rbf_opt_gens, num_interpolants, num_infill_points,
            trace, infill_funcs, infill_iterations, categorical = options

    if trace
        println("Finding new infill samples ...")
    end
    
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
                    MaxFuncEvals=infill_iterations,
                    MaxStepsWithoutProgress=20_000,TraceMode=:silent); 
                    
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
        bo1_solution = BlackBoxOptim.params(pf[idx_obj1]) # get the solution candidate itself... 

        v = copy(permutedims(bo1_solution'))        
        if categorical
            for (i,x_val) in enumerate(v)
                c_ind = closest_index(x_val, vals[i])
                v[i] = vals[i][c_ind]
            end
        end

        # Add the infill point if it does not exist in the plan or infill_plan
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,infill_funcs[i])
        end
    end

    # Add additional points to fill up the desired number of infill points   
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
        if categorical
            for (i,x_val) in enumerate(v)
                c_ind = closest_index(x_val, vals[i])
                v[i] = vals[i][c_ind]
            end
        end

        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,:pareto)
        end
    end
    
    # Add random points if there aren't enough in the pareto front
    while (num_infill_points-size(infill_plan,2)) != 0

        # Add the infill point if it does not exist in the plan or infill_plan
        v = Array{Float64,2}(undef,size(infill_plan,1),1)
        for i = 1:length(v)
            v[i] = _scale(rand(),sr[i][1],sr[i][2],old_min=0,old_max=1)
        end

        # Add the infill point if it does not exist in the plan or infill_plan 
        if categorical
            for (i,x_val) in enumerate(v)
                c_ind = closest_index(x_val, vals[i])
                v[i] = vals[i][c_ind]
            end
        end
        if !colinmat(infill_plan,v) && !colinmat(plan,v)
            push!(infill_prediction,median(sm_interpolant(vec(v))))
            infill_plan = [infill_plan v]            
            push!(infill_type,:rand)
        end
    end

    if trace
        println("Infill samples:")
        println("---------------------------------------------------------------")
        for i = (length(infill_type)-num_infill_points+1):length(infill_type)
            print(@sprintf("%-15s",infill_type[i]))
        end
        print("\n")
        println("---------------------------------------------------------------")
        for j = 1:size(infill_plan,1)
            for i = 1:size(infill_plan,2)
                print(@sprintf("%-15.7g",infill_plan[j,i]))
            end
            print("\n")
        end
        println("---------------------------------------------------------------")
        _, min_loc = findmin(infill_prediction)
        for i = 1:size(infill_plan,2)
            if i == min_loc
                printstyled(@sprintf("%-15.7g",infill_prediction[i]); color=:light_green, bold=true)
            else
                printstyled(@sprintf("%-15.7g",infill_prediction[i]))
            end
        end
        print("\t prediction\n")
        println("---------------------------------------------------------------")
    end

    return infill_plan, criteria, infill_type, infill_prediction
end












function model_infill_brute(sr,plan,samples,sm_interpolant,vals,criteria,options)

    @unpack rbf_opt_gens, num_interpolants, num_infill_points,
            trace, infill_funcs, infill_iterations, categorical = options

    if trace
        println("Finding new infill samples ...")
    end
    
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
    
    num_infill_obj_funs = length(infill_obj_fun(plan[:,1]))
    

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
        fun = functions_to_call[i]
        possible_designs = vec(collect(Base.Iterators.product(vals...)))        

        fun_out = [fun([design...]) for design in possible_designs]        
        
        min_locs = sortperm(fun_out)        
        found_min = false
        j = 1
        while found_min == false            
            v = [possible_designs[min_locs[j]]...]

            # Add the infill point if it does not exist in the plan or infill_plan
            if !colinmat(infill_plan,v) && !colinmat(plan,v)
                push!(infill_prediction,median(sm_interpolant(vec(v))))
                infill_plan = [infill_plan v]            
                push!(infill_type,infill_funcs[i])
                found_min=true
            end
            j+=1
        end
    end

    if trace
        println("Infill samples:")
        println("---------------------------------------------------------------")
        for i = (length(infill_type)-num_infill_points+1):length(infill_type)
            print(@sprintf("%-15s",infill_type[i]))
        end
        print("\n")
        println("---------------------------------------------------------------")
        for j = 1:size(infill_plan,1)
            for i = 1:size(infill_plan,2)
                print(@sprintf("%-15.7g",infill_plan[j,i]))
            end
            print("\n")
        end
        println("---------------------------------------------------------------")
        _, min_loc = findmin(infill_prediction)
        for i = 1:size(infill_plan,2)
            if i == min_loc
                printstyled(@sprintf("%-15.7g",infill_prediction[i]); color=:light_green, bold=true)
            else
                printstyled(@sprintf("%-15.7g",infill_prediction[i]))
            end
        end
        print("\t prediction\n")
        println("---------------------------------------------------------------")
    end

    return infill_plan, criteria, infill_type, infill_prediction
end
