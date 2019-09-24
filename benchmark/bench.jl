# Script used to benchmark the 18 test functions in the cited paper.
# Runs the optimisation for each function 100 times.

using Distributed
addprocs(Sys.CPU_THREADS÷2)   # Set to the number of physical cores in your system (assumes multithreaded CPU)
@everywhere using Statistics
@everywhere using SurrogateModelOptim
@everywhere using ScatteredInterpolation
@everywhere using DelimitedFiles
@everywhere using Distances
@everywhere include(joinpath(dirname(pathof(SurrogateModelOptim)), "../examples/test_functions.jl"))

@everywhere function save_run_bench(fun,test_fun_name,noise,noise_name,path_name)
    func = fun.fun
    search_range = fun.sr
    fun_max = fun.max_val
    fun_min = fun.min_val
    noiselevel = (fun_max-fun_min)*noise

    noisy_fun = function (x)
        func(x)+noiselevel*rand()
    end
    opts = SurrogateModelOptim.Options(SurrogateModelOptim.paper_bench_opts(); parallel_surrogate=false)

    function f()
        result = smoptimize(noisy_fun, search_range;
                        options=opts);

        accumulate(min,[result.infill_samples result.lhc_samples];dims=2)
    end

    res = pmap((x)->f(),1:100)
    res = reduce(vcat,res)

    #Scale results to 0 and 1
    scaled_res = SurrogateModelOptim.scale(res,2,0.0,1.0;old_min=fun_min,old_max=(fun_max+noiselevel))
    
    #Save the results in current folder
    writedlm(joinpath(path_name,test_fun_name*noise_name*".csv"),scaled_res)
end


test_fun_names=[
    "rosenbrock_2D", "rastrigin_2D", "permdbeta_2D", "beale_2D", "goldsteinPrice_2D",
    "sphereFunction_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D",
    "rosenbrock_2D", "rastrigin_2D", "permdbeta_2D", "beale_2D", "goldsteinPrice_2D",
    "sphereFunction_2D", "styblinskiTang_2D", "hart_6D", "rosenbrock_12D"
]

noise = zeros(length(test_fun_names))
noise[length(test_fun_names)÷2+1:end] .= 0.1
noise_name=[["_no_noise" for _ = 1:9]; ["_noise" for _ = 1:9]]
path_name = dirname(@__FILE__)

for i = 1:length(test_fun_names)
    fun = test_funs[Symbol(test_fun_names[i])]
    test_fun_name = test_fun_names[i]
    noise_lvl = noise[i]
    noise_n = noise_name[i]

    save_run_bench(fun,test_fun_name,noise_lvl,noise_n,path_name)
end