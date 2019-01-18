using HypothesisTests





function compare_de_to_surrogate()
    f = rosenbrock2d

    # res = smoptimize(f, [(-5.0,5.0),(-5.0,5.0)],
    #     SurrogateModelOptim.options(num_interpolants=20,
    #                                 num_start_samples=10,
    #                                 sampling_plan_opt_gens=100_000,
    #                                 rbf_opt_gens=30_000,
    #                                 infill_iterations = 25_000,
    #                                 variable_kernel_width = true,
    #                                 variable_dim_scaling = true,
    #                                 num_infill_points = 4,
    #                                 iterations = 1,
    #                                 max_smooth = 100.0,
    #                                 infill_funcs = [:median,:med_2std,:std,:dist],
    #                                 smooth = :single,
    #                                 ));display(funcplotX(x->median(res1.sm_interpolant(x)),-5,5,-5,5))

    

    sol_hist_fitness = Array{Float64,1}()
    sol_hist_iteration = Array{Int64,1}()
    f_c = function (oc)
        push!(sol_hist_fitness, best_fitness(oc))
        push!(sol_hist_iteration, BlackBoxOptim.num_func_evals(oc))
    end        

    res = bboptimize(f; SearchRange=sr,
                        PopulationSize=20, MaxFuncEvals=10000,
                        CallbackFunction = f_c,CallbackInterval = 0.0, 
                        TraceMode = :silent);
    
    equal_iterations_ind = findfirst(x-> x >= length(samples),sol_hist_iteration)
    equal_iterations_fitness = sol_hist_fitness[equal_iterations_ind]

    iterations_ind_for_equal_performance = findfirst(x-> x <= minimum(samples),sol_hist_fitness)
    if iterations_ind_for_equal_performance == nothing
        iterations_for_equal_performance = Inf
    else    
        iterations_for_equal_performance = sol_hist_iteration[iterations_ind_for_equal_performance]
    end

    # println("Function evaluations, ", length(samples), ".  Best DE = ", repr(equal_iterations_fitness), ". Best surrogate = ", minimum(samples))
    # println("DE iterations for >= performance = ", repr(iterations_for_equal_performance), ".")
    # println("Performance increase = ", repr(round(iterations_for_equal_performance/length(samples); digits=2)), ".")
end










# #v1.0 compatible


# using ScatteredInterpolation
# using PGFPlotsX
# using Statistics


# cd("/home/urquhart/Work/Research/SecondProject/PODOptimMethod/podmethodpaper")
# if !in("src/",LOAD_PATH)
#     push!(LOAD_PATH, "src/")
# end
# using RBFoptim_v1



# function fitForrester()
#     #kerns = [Gaussian, Multiquadratic, InverseQuadratic, InverseMultiquadratic]
#     kerns = [Gaussian, InverseQuadratic, InverseMultiquadratic]

#     Forrester(x) = (6*x - 2)^2 * sin(12*x - 4)

#     x = permutedims(collect(range(0, length=5,stop=1)))
#     xLong = permutedims(collect(range(0, length=100,stop=1)))
#     y = Forrester.(x)
#     yLong = Forrester.(xLong)
#     runs = 1000 #1000 in paper
#     out = Array{Float64}(undef,runs,length(xLong))
#     i = 1
#     while i <= runs

        
#         result = RBF_hypers_opt(y,x,kerns,10000;
#             rippa = true,
#             fixedKernel_and_Width = false,
#             fixedDimScaling = true,
#             TraceMode = :none,
#             condmax = 1e4,
#             maxWidth = 1000.0,
#             maxScale = 10.0)

#         if true#result[1].fitness <= 1e-10
#             out[i,:] = evalInterpolant(
#                 result[1],
#                 x,
#                 y',
#                 xLong,
#                 minimum(x,dims=2),
#                 maximum(x,dims=2)
#                 )
#             @show i+=1    
#         end
#     end
#     x,y,out,xLong,yLong
# end

# x,y,out,xLong,yLong=fitForrester()


# #ymeanplot = vec(mean(out,dims=1))
# ymeanplot = vec(median(out,dims=1))
# xplot = vec(xLong)
# ystdplot = vec(std(out,dims=1))


# data = rand(Normal(10,10),100)
# ci(OneSampleTTest(data))

# if !in("\\usepgfplotslibrary{fillbetween}",PGFPlotsX.CUSTOM_PREAMBLE)
#     push!(PGFPlotsX.CUSTOM_PREAMBLE,"\\usepgfplotslibrary{fillbetween}")
# end

# ax = @pgf Axis(
#     {
#         #xmajorgrids,
#         #ymajorgrids,
#         "axis x line=bottom",
#         "axis y line=left",
#         "legend style" =
#             {
#                 "draw=none",
#                 "at={(0.03,.97)}",
#                 "anchor=north west",
#                 "font=\\footnotesize",
#             },
#     },
#     Plot(
#         {
#             no_marks,
#             color = "black",
#         },
#         Coordinates(xplot, vec(yLong))
#     ),
#     LegendEntry("Forrester function"),
#     PlotInc(
#             {
#                 "only marks",
#                 mark = "*",
#                 "mark size" = 1.7,
#                 color = "black",
#                 "mark options" = {fill="black"},
#             },
#             Coordinates(vec(x),vec(y))
#         ),
#     LegendEntry("Interpolation points"),
#     PlotInc(
#         {
#             no_marks,
#             color = "black",
#             dashed,
#         },
#         Coordinates(xplot, ymeanplot)
#     ),
#     LegendEntry("Average interpolant"),
#     PlotInc(
#         {
#             "name path=upper",
#             "draw=none",
#             no_marks,
#             "forget plot",
#         },
#         Table([:x => xplot, :y => ymeanplot.+ystdplot]),
#     ),
#     PlotInc(
#         {
#             "name path=lower",
#             "draw=none",
#             no_marks,
#             "forget plot",
#         },
#         Table([:x => xplot, :y => ymeanplot.-ystdplot]),
#     ),
#     PlotInc(
#         {
#             "name path=upper2",
#             "draw=none",
#             no_marks,
#             "forget plot",
#         },
#         Table([:x => xplot, :y => ymeanplot.+2*ystdplot]),
#     ),
#     PlotInc(
#         {
#             "name path=lower2",
#             "draw=none",
#             no_marks,
#             "forget plot",
#         },
#         Table([:x => xplot, :y => ymeanplot.-2*ystdplot]),
#     ),
    
#     "\\addplot [fill=black!20, area legend] fill between[of=upper2 and lower2];",
#     LegendEntry("\$\\pm 2\$ Standard deviation"),
#     "\\addplot [fill=black!40, area legend] fill between[of=upper and lower];",
#     LegendEntry("\$\\pm 1\$ Standard deviation"),
# )

# #pgfsave("paper/figures/elephant.tikz", ax;dpi=300)






# using ColorSchemes
# cols = ColorSchemes.Blues_9;
# ax = @pgf Axis(
#     {
#         #xmajorgrids,
#         #ymajorgrids,
#         "axis x line=bottom",
#         "axis y line=left",
#         "legend style" =
#             {
#                 "draw=none",
#                 "at={(0.03,.97)}",
#                 "anchor=north west",
#                 "font=\\footnotesize",
#             },
#     },
#     PlotInc(
#             {
#                 "only marks",
#                 mark = "*",
#                 "mark size" = 1.7,
#                 color = "black",
#                 "mark options" = {fill="black"},
#             },
#             Coordinates(vec(x),vec(y))
#         ),
#     LegendEntry("Interpolation points"),
#     PlotInc(
#             {
#                 no_marks,
#                 color = "black!50!white",
#                 opacity = "0.4",
#                 "line width=0.5pt",
#             },
#     Coordinates(xplot, out[1,:]),
#     ),
#     LegendEntry("Surrogate "),
#     PlotInc(
#         {
#             no_marks,
#             color = cols[8],
#             "line width=1pt",
#         },
#         Coordinates(xplot, vec(yLong))
#     ),
#     LegendEntry("Forrester function"),
# );

# for i = 2:size(out,1)
#     p = @pgf Plot(
#         {
#             no_marks,
#             color = "black!50!white",
#             solid,
#             opacity = "0.4",
#             "line width=0.5pt",
#         },
#         Coordinates(xplot, out[i,:]),
#     )
#     push!(ax,p)
# end
# p = @pgf PlotInc(
#     {
#         no_marks,
#         color = cols[8],
#         "line width=1pt",
#     },
#     Coordinates(xplot, vec(yLong))
# );
# push!(ax,p);
# ax

# #pgfsave("paper/figures/elephantlines.tikz", ax;dpi=300)
