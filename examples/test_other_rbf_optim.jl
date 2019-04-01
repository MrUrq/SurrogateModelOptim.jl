using Distributed
using DelimitedFiles
using Parameters
using Clustering
using SurrogateModelOptim
using Distances
using Optim
using BlackBoxOptim
using Statistics

options = SurrogateModelOptim.Options(
    num_interpolants=20, #Preferably even number of added processes
    rbf_opt_gens=50_000, infill_iterations=25_000,
    num_infill_points=4, trace=true,
    infill_funcs = [:std,:mean,:median,:min],
    rbf_opt_method = :de_rand_1_bin)


# Read in the plan for the successful runs with their corresponding Cl/Cl
file_name = "/home/urquhart/Work/Courses_teaching/2019_RVAD/DiffuserOptim/V2/successful_runs.csv"
plan_vals = readdlm(file_name,',')[2:end,4:end-4]
samples_org = permutedims(readdlm(file_name,',')[2:end,3].*-1)[:,1:100]   

#"D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"
plan = Array{Float64,2}(permutedims(plan_vals))[:,1:100]






@unpack rippa, variable_kernel_width, variable_dim_scaling, rbf_opt_method, 
            min_rbf_width, max_rbf_width, min_scale, max_scale, cond_max,
            rbf_dist_metric, rbf_opt_gens, kerns, smooth, max_smooth, smooth_user = options

samples = vec(samples_org)


# Create the hyperparameter search range based on the input options 
sr = SurrogateModelOptim.construct_search_range(plan, variable_kernel_width,
                            min_rbf_width, max_rbf_width, variable_dim_scaling,
                            min_scale, max_scale, smooth, max_smooth)

# RBF hyperparameter objective function
itp_obj = function (x)
    x_max = [x for (_,x) in sr]
    x_min = [x for (x,_) in sr]

    for (i,xv) in enumerate(x)
        (xv > x_max[i]) && (x[i] = x_max[i])
        (xv < x_min[i]) && (x[i] = x_min[i])
    end
    SurrogateModelOptim.interp_obj(x,kerns,samples,plan; 
            rippa=rippa, variable_kernel_width=variable_kernel_width,
            variable_dim_scaling=variable_dim_scaling, smooth=smooth,
            cond_max=cond_max)
end



# Optimize the interpolant hyperparameters
function de(PopulationSize)
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=rbf_opt_gens,
            TraceMode=:verbose, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6,
            PopulationSize = PopulationSize,
            MaxStepsWithoutProgress=rbf_opt_gens,
            MaxNumStepsWithoutFuncEvals=rbf_opt_gens,
            )
    return res.archive_output.best_fitness
end

function ne()
    x0 = [rand()*(y-x)+x for (x,y) in sr]
    r = optimize(itp_obj,x0,NelderMead(),
            Optim.Options(  f_calls_limit = rbf_opt_gens,
                            iterations = rbf_opt_gens,                            
                            g_tol = 1e-12,
                            show_trace = true,
                            allow_f_increases = true))
    return r.minimum
end

function lbfgs()
    x0 = [rand()*(y-x)+x for (x,y) in sr]
    r = optimize(itp_obj,x0,LBFGS(),
            Optim.Options(  f_calls_limit = rbf_opt_gens,
                            g_calls_limit = g_calls_limit,
                            iterations = rbf_opt_gens,
                            g_tol = 1e-12,
                            show_trace = true,
                            allow_f_increases = true))
    return r.minimum
end

function ne_hybrid()
    per = 0.2
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=round(Int,per*rbf_opt_gens),
            TraceMode=:verbose, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6);

    x0 = best_candidate(res)
    r = optimize(itp_obj,x0,NelderMead(),
            Optim.Options(  f_calls_limit = round(Int,(1-per)*rbf_opt_gens),
                            iterations = rbf_opt_gens,
                            g_tol = 1e-12,
                            show_trace = true,
                            allow_f_increases = true))
    return r.minimum
end

function lbfgs_hybrid()
    per = 0.75
    res = bboptimize(itp_obj; 
            Method=rbf_opt_method,SearchRange=sr, MaxFuncEvals=round(Int,per*rbf_opt_gens),
            TraceMode=:verbose, rbf_dist_metric=rbf_dist_metric,
            TargetFitness = 1e-5, FitnessTolerance = 1e-6);

    x0 = best_candidate(res)
    r = optimize(itp_obj,x0,ParticleSwarm(),
            Optim.Options(  f_calls_limit = round(Int,(1-per)*rbf_opt_gens),
                            iterations = rbf_opt_gens,
                            g_tol = 1e-12,
                            show_trace = true,
                            allow_f_increases = true))
    return r.minimum
end

t1 = map(x->de(5),1:20)
t2 = map(x->de(25),1:20)  #winner for 25, 50 and 100 sample points
t3 = map(x->de(50),1:20)
t4 = map(x->de(75),1:20)
t5 = map(x->de(100),1:20)
#t2 = map(x->ne(),1:20)
#t3 = map(x->hybrid(),1:20)

#t4 = map(x->lbfgs_hybrid(),1:1)