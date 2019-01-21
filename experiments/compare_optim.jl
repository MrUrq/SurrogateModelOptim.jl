using Distributions

function conf(conf_level,d)
    l = length(d)
    f_mean = mean(d)

    α = (1 - conf_level)
    tstar = quantile(TDist(l-1), 1 - α/2)
    SE = std(d; mean = f_mean)/sqrt(l)
    c = tstar * SE
    return (f_mean, f_mean - c, f_mean + c) 
end



function conf_subtract(conf_level,d1,d2)

    l1 = length(d1)
    l2 = length(d2)
    f_mean1 = mean(d1)
    f_mean2 = mean(d2)
    
    α = (1 - conf_level)
    tstar = quantile(TDist(l-1), 1 - α/2)
    std1 = std(d1; mean = f_mean1)
    std2 = std(d2; mean = f_mean2)

    f_m = f_mean1-f_mean2
    c = tstar * sqrt(std1^2/l1 + std2^2/l2)
    return (f_m , f_m - c, f_m + c) 
end

function conf_division(conf_level,d1,d2)

    l1 = length(d1)
    l2 = length(d2)
    f_mean1 = mean(d1)
    f_mean2 = mean(d2)
    
    α = (1 - conf_level)
    tstar = quantile(TDist(l-1), 1 - α/2)
    std1 = std(d1; mean = f_mean1)
    std2 = std(d2; mean = f_mean2)

    f_m = f_mean1/f_mean2
    c = tstar * abs(f_m) * sqrt(std1^2/f_mean1/l1 + std2^2/f_mean2/l2)
    return (f_m , f_m - c, f_m + c) 
end

function optim_res_compare_conf(conf_level,swept_runs1,swept_runs2)
    
    x1, mean_opt1, _, _ = optim_res_to_conf(conf_level,swept_runs1)
    x2, mean_opt2, _, _ = optim_res_to_conf(conf_level,swept_runs2)
    
    comparison_x = x1[[any(x .== x2) for x in x1]]

    f_evals = length(comparison_x)
    mean_opt = Array{Float64,1}(undef,f_evals)
    mean_opt_lo = Array{Float64,1}(undef,f_evals)
    mean_opt_hi = Array{Float64,1}(undef,f_evals)

    for (i,ind) in enumerate(comparison_x)
        
        d1 = collect(skipmissing(swept_runs1[:,ind]))
        d2 = collect(skipmissing(swept_runs2[:,ind]))
        
        mean_opt[i], mean_opt_lo[i], mean_opt_hi[i] = 1 .- conf_division(conf_level,d1,d2)
    end

    return comparison_x, mean_opt, mean_opt_lo, mean_opt_hi
end

function optim_res_to_conf(conf_level,swept_runs)
    f_evals = size(swept_runs,2)
    mean_opt = Array{Union{Missing, Float64},1}(undef,f_evals)
    mean_opt_lo = Array{Union{Missing, Float64},1}(undef,f_evals)
    mean_opt_hi = Array{Union{Missing, Float64},1}(undef,f_evals)

    for i = 1:size(swept_runs,2)

        data = swept_runs[:,i]        
        d = collect(skipmissing(data))

        if length(d) < 2
            mean_opt[i] = missing
            mean_opt_lo[i] = missing
            mean_opt_hi[i] = missing
        else
            mean_opt[i], mean_opt_lo[i], mean_opt_hi[i] = conf(conf_level,d)
        end
    end

    x = findall(!ismissing,mean_opt)
    mean_opt = Array{Float64}(mean_opt[x])
    mean_opt_lo = Array{Float64}(mean_opt_lo[x])
    mean_opt_hi = Array{Float64}(mean_opt_hi[x])

    return x, mean_opt, mean_opt_lo, mean_opt_hi
end





####################################### Walk backwards adding the last seen value
##################### Push every function evaluation to a new array to take it out from DE optim.
function de_stats(f)

    num_runs = 1000
    f_evals = 10000
    sr = [(-5.0,5.0),(-5.0,5.0)]
    swept_runs = Array{Union{Missing, Float64},2}(undef,num_runs,f_evals).=missing
    for i = 1:num_runs
        
        f_c = function (oc)
            swept_runs[i,BlackBoxOptim.num_func_evals(oc)] = best_fitness(oc)
        end        

        res = bboptimize(f; SearchRange=sr,
                            PopulationSize=20, MaxFuncEvals=f_evals,
                            CallbackFunction = f_c,CallbackInterval = 0.0, 
                            TraceMode = :silent);
    end

    swept_runs = SurrogateModelOptim._scale.(swept_runs,0,1,old_min=min_val,old_max=max_val)

    x, mean_opt, mean_opt_lo2, mean_opt_hi2 = optim_res_to_conf(0.66,swept_runs)
    _, mean_opt, mean_opt_lo, mean_opt_hi = optim_res_to_conf(0.95,swept_runs)

    return swept_runs, x, mean_opt, mean_opt_lo2, mean_opt_hi2, mean_opt_lo, mean_opt_hi
end

using JLD
function sm_stats(test_f,num_runs,f_evals,num_lhc_points)

    sr = test_f.sr
    f = test_f.fun
    min_val = test_f.min_val
    max_val = test_f.max_val

    iterations = f_evals-num_lhc_points

    swept_runs = Array{Union{Missing, Float64},2}(undef,num_runs,f_evals).=missing

    for i = 1:num_runs

        @show i

        res = smoptimize(f, sr,
            SurrogateModelOptim.options(num_interpolants=20,
                                        num_start_samples=num_lhc_points,
                                        sampling_plan_opt_gens=100_000,
                                        rbf_opt_gens=25_000,
                                        infill_iterations = 25_000,
                                        variable_kernel_width = true,
                                        variable_dim_scaling = true,
                                        num_infill_points = 1,
                                        iterations = iterations,
                                        max_smooth = 100.0,
                                        infill_funcs = [:median,:med_2std],
                                        smooth = :single,
                                        ))
        
        out = Array{Float64,2}(undef,1,f_evals)
        samples = SurrogateModelOptim._scale(vec([res.lhc_samples res.infill_samples]),0,1,old_min=min_val,old_max=max_val)
        tmp_sample = Inf
        for i = 1:length(out) 
            if samples[i] < tmp_sample
                out[i] = samples[i]
                tmp_sample = samples[i]
            else
                out[i] = out[i-1]
            end     
        end
        
        swept_runs[i,:] = out
    end
    
    x, mean_opt, mean_opt_lo66, mean_opt_hi66 = optim_res_to_conf(0.66,swept_runs)
    _, _, mean_opt_lo95, mean_opt_hi95 = optim_res_to_conf(0.95,swept_runs)

    
    save("/home/urquhart/.julia/dev/SurrogateModelOptim/experiments/my_test_file.jld", "swept_runs", swept_runs)
    
    return swept_runs, x, mean_opt, mean_opt_lo66, mean_opt_hi66, mean_opt_lo95, mean_opt_hi95
end




# mean_opt2 = mean_opt[mean_opt .> 0.001]
# x = x[mean_opt .> 0.001]
# mean_opt_lo = mean_opt_lo[mean_opt .> 0.001]
# mean_opt2 = mean_opt[mean_opt .> 0.001]
# x = x[mean_opt .> 0.001]
# mean_opt_lo = mean_opt_lo[mean_opt .> 0.001]
# mean_opt_hi = mean_opt_hi[mean_opt .> 0.001]
# mean_opt_lo2 = mean_opt_lo2[mean_opt .> 0.001]
# mean_opt_hi2 = mean_opt_hi2[mean_opt .> 0.001]
# mean_opt = mean_opt2



# mean_opt_hi = mean_opt_hi[mean_opt .> 0.001]
# mean_opt2 = mean_opt[mean_opt .> 0.001]
# x = x[mean_opt .> 0.001]
# mean_opt_lo = mean_opt_lo[mean_opt .> 0.001]
# mean_opt_hi = mean_opt_hi[mean_opt .> 0.001]
# mean_opt_lo2 = mean_opt_lo2[mean_opt .> 0.001]
# mean_opt_hi2 = mean_opt_hi2[mean_opt .> 0.001]
# mean_opt = mean_opt2



# mean_opt_lo2 = mean_opt_lo2[mean_opt .> 0.001]
# mean_opt_hi2 = mean_opt_hi2[mean_opt .> 0.001]
# mean_opt = mean_opt2






using PGFPlotsX
if !in("\\usepgfplotslibrary{fillbetween}",PGFPlotsX.CUSTOM_PREAMBLE)
    push!(PGFPlotsX.CUSTOM_PREAMBLE,"\\usepgfplotslibrary{fillbetween}")
end

function plot_confint(x,mean_opt_org,mean_opt_hi66_org,mean_opt_lo66_org,mean_opt_hi95_org,mean_opt_lo95_org,num_lhc_points)

    mean_opt_hi66 = copy(mean_opt_hi66_org)
    mean_opt_lo66 = copy(mean_opt_lo66_org)
    mean_opt_hi95 = copy(mean_opt_hi95_org)
    mean_opt_lo95 = copy(mean_opt_lo95_org)

    min_val = minimum([ mean_opt_hi66[mean_opt_hi66 .> 0]
                            mean_opt_lo66[mean_opt_lo66 .> 0]
                            mean_opt_hi95[mean_opt_hi95 .> 0]
                            mean_opt_lo95[mean_opt_lo95 .> 0]])
    
    mean_opt_hi66[mean_opt_hi66 .< 0] .= min_val
    mean_opt_lo66[mean_opt_lo66 .< 0] .= min_val
    mean_opt_hi95[mean_opt_hi95 .< 0] .= min_val
    mean_opt_lo95[mean_opt_lo95 .< 0] .= min_val

    
    ax = @pgf Axis(
        {
            #xmajorgrids,
            #ymajorgrids,
            "axis x line=bottom",
            "axis y line=left",
            "legend style" =
                {
                    "draw=none",
                    "at={(0.99,.97)}",
                    "anchor=north east",
                    "font=\\footnotesize",
                },
            ymode = "log"
        },
        PlotInc(
            {
                no_marks,
                color = "black",
                dashed,
            },
            Coordinates(x, mean_opt)
        ),
        LegendEntry("Average interpolant"),
        PlotInc(
            {
                "name path=upper66",
                "draw=none",
                no_marks,
                "forget plot",
            },
            Table([:x => x, :y => mean_opt_hi66]),
        ),
        PlotInc(
            {
                "name path=lower66",
                "draw=none",
                no_marks,
                "forget plot",
            },
            Table([:x => x, :y => mean_opt_lo66]),
        ),
        PlotInc(
            {
                "name path=upper95",
                "draw=none",
                no_marks,
                "forget plot",
            },
            Table([:x => x, :y => mean_opt_hi95]),
        ),
        PlotInc(
            {
                "name path=lower95",
                "draw=none",
                no_marks,
                "forget plot",
            },
            Table([:x => x, :y => mean_opt_lo95]),
        ),
        "\\addplot [fill=black!20, area legend] fill between[of=upper95 and lower95];",
        LegendEntry("\$\\pm 2\$ Standard deviation"),
        "\\addplot [fill=black!40, area legend] fill between[of=upper66 and lower66];",
        LegendEntry("\$\\pm 1\$ Standard deviation"),
        PlotInc(
            {
                no_marks,
                color = "black",
                solid,
                "forget plot",
            },
            Coordinates([10, 10], [1, min_val])
        ),
    )
end

#pgfsave("paper/figures/elephant.tikz", ax;dpi=300)

x,mean_opt,mean_opt_lo66,mean_opt_hi66=optim_res_compare_conf(0.66,swept_runs,swept_runs_de)
x,mean_opt,mean_opt_lo95,mean_opt_hi95=optim_res_compare_conf(0.95,swept_runs,swept_runs_de)


x, mean_opt, mean_opt_lo66, mean_opt_hi66 = optim_res_to_conf(0.66,swept_runs)
       _, _, mean_opt_lo95, mean_opt_hi95 = optim_res_to_conf(0.95,swept_runs)
plot_confint(x,mean_opt,mean_opt_hi66,mean_opt_lo66,mean_opt_hi95,mean_opt_lo95,num_lhc_points)

min_val = minimum(mean_opt)

x, mean_opt, mean_opt_lo66, mean_opt_hi66 = optim_res_to_conf(0.66,swept_runs_de[:,1:25])
       _, _, mean_opt_lo95, mean_opt_hi95 = optim_res_to_conf(0.95,swept_runs_de[:,1:25])
plot_confint(x,mean_opt,mean_opt_hi66,mean_opt_lo66,mean_opt_hi95,mean_opt_lo95,num_lhc_points)



using ColorSchemes

function plot_all(swept_runs1,swept_runs2)
    cols = ColorSchemes.Blues_9;

    x1, mean_opt1, _, _ = optim_res_to_conf(0.95,swept_runs1)
    x2, mean_opt2, _, _ = optim_res_to_conf(0.95,swept_runs2)
    
    x = x1[[any(x .== x2) for x in x1]]

    med1 = [median(collect(skipmissing(swept_runs1[:,i]))) for i in x]
    med2 = [median(collect(skipmissing(swept_runs2[1:10,i]))) for i in x]
    

    ax = @pgf Axis(
        {
            #xmajorgrids,
            #ymajorgrids,
            "axis x line=bottom",
            "axis y line=left",
            "legend style" =
                {
                    "draw=none",
                    "at={(0.03,.97)}",
                    "anchor=north west",
                    "font=\\footnotesize",
                },
            ymode = "log",
        },
        PlotInc(
                {
                    no_marks,
                    color = "black",
                    "line width=1pt",
                },
        Coordinates(x, med1),
        ),
        LegendEntry("Median 1"),

        PlotInc(
                {
                    no_marks,
                    color = "black",
                    "line width=1pt",
                },
        Coordinates(x, med2),
        ),
        LegendEntry("Median 2"),
    );

    
    for i = 1:10
        xp = Array{Int64,1}()
        yp = Array{Float64,1}()

        for j = 1:25
            val = swept_runs2[i,j]
            if !ismissing(val)
                push!(xp,j) 
                push!(yp,val)
            end
        end

        p = @pgf PGFPlotsX.Plot(
            {
                no_marks,
                color = "red",
                solid,
                opacity = "0.4",
                "line width=0.5pt",
            },
            Coordinates(xp, yp),
        )
        push!(ax,p)
    end


    for i = 1:10
        xp = Array{Int64,1}()
        yp = Array{Float64,1}()

        for j = 1:25
            val = swept_runs1[i,j]
            if !ismissing(val)
                push!(xp,j) 
                push!(yp,val)
            end
        end

        p = @pgf PGFPlotsX.Plot(
            {
                no_marks,
                color = "blue",
                solid,
                opacity = "0.4",
                "line width=0.5pt",
            },
            Coordinates(xp, yp),
        )
        push!(ax,p)
    end

    ax
end
plot_all(swept_runs,swept_runs_de)
# #pgfsave("paper/figures/elephantlines.tikz", ax;dpi=300)
