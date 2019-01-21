using PGFPlotsX
using Colors
using ColorSchemes
using Random
using SurrogateModelOptim
using ColorBrewer
using Distances
using ScatteredInterpolation
using Statistics

function funcplotX(fun,xmin,xmax,ymin,ymax)
    
    np = 40   #Number of points in 1 dimension for plotting
    x = range(xmin,stop=xmax,length=np)
    y = range(ymin,stop=ymax,length=np)

    asd = (x,y)->[x,y]
    data = fun.(asd.(x,y'))
    
    ax = @pgf Axis(
        {
            #view = (-24, 40),
            #"point meta min" = (minimum(data)),
            #"point meta max" = (maximum(data)*0.75),
            #zmin = -20,
            #zmax = 2020,
        },
        Plot3(
            {
                surf,
            },
            Coordinates(x,y,data)
        )
    )
    @pgf push!(ax, ("PuBuGn_9", sortcolorscheme(ColorSchemes.PuBuGn_9, rev=false)))    

    return ax
end



test_f=test_funs[:rosenbrock_2D]

sr = test_f.sr
f = test_f.fun

res1 = smoptimize(f, sr,	
        SurrogateModelOptim.options(num_interpolants=20,
                                    num_start_samples=10,	           
                                    sampling_plan_opt_gens=100_000,	   
                                    rbf_opt_gens=30_000,	           
                                    infill_iterations = 25_000,	       
                                    variable_kernel_width = true,	   
                                    variable_dim_scaling = true,	   
                                    num_infill_points = 4,	           
                                    iterations = 15,	               
                                    #max_smooth = 100.0,	               
                                    infill_funcs = [:median,:med_2std,:std,:dist],	         
                                    smooth = :single,	                                   
                                    ));display(funcplotX(x->median(res1.sm_interpolant(x)),-5,5,-5,5))
