using PGFPlotsX
using Colors
using ColorSchemes
using Random
using SurrogateModelOptim
using ColorBrewer
using Distances
using ScatteredInterpolation
using Statistics

#ax = funcplotX(x->median(asd[3](x))[1],-5,5,-5,5)

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


ax = funcplotX(rosenbrock2d,-5,5,-5,5)

# asd = smoptimize(x->(styblinskiTang(x)+0*rand(MersenneTwister(abs(sum(reinterpret(Int64,x)))))), [(-5.0,5.0),(-5.0,5.0)],SurrogateModelOptim.options(num_interpolants=20, num_start_samples=20, sampling_plan_opt_gens=100000, rbf_opt_gens=10000,  variable_kernel_width = true,
#                                                                                                 variable_dim_scaling = true,
#                                                                                                 smooth = false));

@time asd = smoptimize(x->(rosenbrock2d(x)+0*rand(MersenneTwister(abs(sum(reinterpret(Int64,x)))))), [(-5.0,5.0),(-5.0,5.0)],SurrogateModelOptim.options(num_interpolants=20, num_start_samples=10, sampling_plan_opt_gens=100000, rbf_opt_gens=10000,  variable_kernel_width = true,
                                                                                                       variable_dim_scaling = true,
                                                                                                       ));

ax = funcplotX(x->median(asd[3](x))[1],-5,5,-5,5)

