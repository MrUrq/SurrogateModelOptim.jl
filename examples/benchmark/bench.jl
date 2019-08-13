using Distributed
addprocs(20)
@everywhere using Statistics
@everywhere using SurrogateModelOptim
@everywhere using ScatteredInterpolation
@everywhere using DelimitedFiles
@everywhere using Distances
@everywhere include(raw"/home/urquhart/.julia/dev/SurrogateModelOptim/examples/test_functions.jl")


@everywhere function save_run_bench(fun,test_fun_name,noise,noise_name,smooth,path_name)
    func = fun.fun
    search_range = fun.sr
    fun_max = fun.max_val
    fun_min = fun.min_val
    noiselevel = (fun_max-fun_min)*noise

    noisy_fun = function (x)
        func(x)+noiselevel*rand()
    end

    function f()
        result = smoptimize(noisy_fun, search_range;
                        options=SurrogateModelOptim.Options(
                            num_start_samples = 5,
                            trace = :compact,
                            sampling_plan_opt_gens = 1_000,
                            rippa = true,
                            kerns = [ScatteredInterpolation.Gaussian],
                            rbf_opt_gens = 1_000,
                            rbf_opt_pop = 50,
                            rbf_opt_method = :adaptive_de_rand_1_bin_radiuslimited,
                            rbf_dist_metric = Distances.Euclidean(),
                            variable_kernel_width = true,
                            variable_dim_scaling = true,
                            cond_max = 5e12,
                            cond_check = false,
                            max_rbf_width = 1.0,
                            min_rbf_width = 0.0,
                            max_scale = 1.0,
                            min_scale = 0.0,
                            num_interpolants = 10,
                            smooth = smooth,
                            max_smooth = 0.005,
                            smooth_user = 0.0,
                            iterations = 100,
                            num_infill_points = 1,
                            parallel_surrogate = false,
                            infill_funcs = [:median,:std],
                            infill_iterations = 10_000,
                            create_final_surrogate = false,
                            ));
            
        accumulate(min,[result.infill_samples result.lhc_samples];dims=2)
    end
    res = pmap((x)->f(),1:100)
    res = reduce(vcat,res)

    #Scale results to 0 and 1
    scaled_res = SurrogateModelOptim.scale(res,2,0.0,1.0;old_min=fun_min,old_max=(fun_max+noiselevel))
    
    writedlm(joinpath(path_name,test_fun_name*noise_name*".csv"),scaled_res)
end


test_fun_names=["rastrigin_2D", "rosenbrock_2D", "permdbeta_2D", "beale_2D", "goldsteinPrice_2D", "sphereFunction_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D", "rastrigin_2D", "rosenbrock_2D", "permdbeta_2D", "beale_2D", "goldsteinPrice_2D", "sphereFunction_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D"]

noise = zeros(length(test_fun_names))
noise[length(test_fun_names)÷2+1:end] .= 0.1

noise_name=[["_no_noise" for _ = 1:9]; ["_noise" for _ = 1:9]]

smooth=[[false for _ = 1:9]; [:single for _ = 1:9]]







# test_fun_names=["rastrigin_2D", "rosenbrock_2D", "permdbeta_2D", "beale_2D", "goldsteinPrice_2D", "sphereFunction_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D", "rastrigin_2D", "rosenbrock_2D", "permdbeta_2D", "beale_2D", "goldsteinPrice_2D", "sphereFunction_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D"]

# noise = zeros(length(test_fun_names))
# noise[length(test_fun_names)÷2+1:end] .= 0.1

# noise_name=[["_no_noise" for _ = 1:9]; ["_noise" for _ = 1:9]]

# smooth=[[false for _ = 1:9]; [:single for _ = 1:9]]


test_fun_names = ["styblinskiTang_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D", "hart_6D"]
noise = [0.0,0.1,0.1,0.1,0.0]
noise_name = ["_no_noise","_noise","_noise","_noise","_no_noise_w_smoothing"]
smooth = [false, :single, :single, :single, :single] 

path_name = dirname(@__FILE__)

for i = 1:length(test_fun_names)
    fun = test_funs[Symbol(test_fun_names[i])]
    test_fun_name = test_fun_names[i]
    noise_lvl = noise[i]
    smoothing = smooth[i]
    noise_n = noise_name[i]

    save_run_bench(fun,test_fun_name,noise_lvl,noise_n,smoothing,path_name)
end


