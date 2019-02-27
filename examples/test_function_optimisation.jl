using SurrogateModelOptim
using PlotlyJS
using Statistics


# Styblinski-Tang test function
function styblinskiTang_2D(x)
    out = 0.0    
    xy = [x[1],x[2]]
    for i = 1:2
        out += xy[i]^4 - 16*xy[i]^2 + 5*xy[i]
    end
    out *= 0.5
    return out
end

function rosenbrock_2D(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function rastrigin_ND(x)
    d=length(x)
    sum = 0
    for ii = 1:d
        xi = x[ii]
        sum = sum + (xi^2 - 10*cos(2*pi*xi))
    end
    out = 10*d + sum
end


# Optimize the test function
#sr = [(-5.0,5.0),(-5.0,5.0)]
sr = [(-5.12,5.12) for i = 1:20]
result = smoptimize(rastrigin_ND, sr, SurrogateModelOptim.Options(
                                                iterations=300-20,
                                                num_interpolants=20,
                                                num_start_samples=20,
                                                rbf_opt_gens=50_000,
                                                infill_iterations=100_000,
                                                smooth=:single,
                                                num_infill_points=1,
                                                    ));


# Plot results
function plot_fun_2D(fun,sr)    
    N = 51    
    x = range(sr[1][1], stop = sr[1][2], length = N)
    y = range(sr[2][1], stop = sr[2][2], length = N)

    grid(x,y) = [x,y]
    z = @. fun(grid(x,y'))

    trace = surface(x=x,y=y,z=z, colorscale="Viridis")     
    p = plot(trace) 
end

# Plot error
function plot_error_2D(fun1,fun2,sr)    
    N = 51    
    x = range(sr[1][1], stop = sr[1][2], length = N)
    y = range(sr[2][1], stop = sr[2][2], length = N)

    grid(x,y) = [x,y]
    z1 = @. fun1(grid(x,y'))
    z2 = @. fun2(grid(x,y'))
    z = abs.(z1.-z2)
    
    trace = contour(x=x,y=y,z=z)

    p = plot(trace) 
end

#plot_fun_2D(styblinskiTang_2D,sr)
plot_fun_2D(x->median(result.sm_interpolant(x)),sr)
#plot_error_2D(styblinskiTang_2D,x->median(result.sm_interpolant(x)),sr)

