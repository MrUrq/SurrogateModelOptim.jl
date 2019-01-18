using PlotlyJS
using ColorBrewer

function plot_sm2D(f,fun,sr)

    x, y = range(sr[1][1], stop=sr[1][2], length=100), range(sr[2][1], stop=sr[2][2], length=100)

    Xgrid = [i for i in x, j in y]
    Ygrid = [j for i in x, j in y]

    function plot_fun(fun)
        grid(x,y) = [x,y]
        Y = @. median(fun(grid(Xgrid,Ygrid)))

        col = [[i, reverse(ColorBrewer.colorSchemes["PuBuGn"]["9"])[ind]] for (ind,i) in enumerate(range(0, stop=1, length=7))]

        trace = surface(x=Xgrid,y=Ygrid,z=Y,colorscale=col,reversescale=true) 
        layout = Layout(height=440, 
                        scene = (   xaxis=attr(title="X"), 
                                    yaxis=attr(title="Y"), 
                                    zaxis=attr(title="Z",)),#range=[-10,10])), 
                        margin=attr(l=20, r=20, b=20, t=20), 
                        ) 
        p = plot(trace, layout) 
    end

    [plot_fun(f) plot_fun(fun)]
end

function plot_sm2D(f,sr)

    x, y = range(sr[1][1], stop=sr[1][2], length=100), range(sr[2][1], stop=sr[2][2], length=100)

    Xgrid = [i for i in x, j in y]
    Ygrid = [j for i in x, j in y]

    function plot_fun(fun)
        grid(x,y) = [x,y]
        Y = @. median(fun(grid(Xgrid,Ygrid)))

        col = [[i, reverse(ColorBrewer.colorSchemes["PuBuGn"]["9"])[ind]] for (ind,i) in enumerate(range(0, stop=1, length=7))]

        trace = surface(x=Xgrid,y=Ygrid,z=Y,colorscale=col,reversescale=true) 
        layout = Layout(height=440, 
                        scene = (   xaxis=attr(title="X"), 
                                    yaxis=attr(title="Y"), 
                                    zaxis=attr(title="Z",)),#range=[-10,10])), 
                        margin=attr(l=20, r=20, b=20, t=20), 
                        ) 
        p = plot(trace, layout) 
    end

    plot_fun(f)
end

#plot_sm2D(rosenbrock2d,res1.sm_interpolant,[(-5.0,5.0),(-5.0,5.0)]);



            
