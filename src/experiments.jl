using PGFPlotsX
using Colors
using ColorSchemes
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