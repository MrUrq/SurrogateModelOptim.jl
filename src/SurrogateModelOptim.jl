module SurrogateModelOptim

export  smoptimize,
        rosenbrock2d,
        rotatedHyperElipsoid,
        styblinskiTang,
        hart6,
        minimum,
        maximum,
        std,
        median,
        mean,
        model_infill

function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
function rotatedHyperElipsoid(x)
    out = 0.0
    xy = [x[1],x[2]]
    for i = 1:2
        inner = 0.0
        for j = 1:i
            inner = inner + xy[j]^2
        end
        out += inner
    end
    return out
end
function styblinskiTang(x)
    out = 0.0    
    xy = [x[1],x[2]]
    for i = 1:2
        out += xy[i]^4 - 16*xy[i]^2 + 5*xy[i]
    end
    out *= 0.5
    return out
end



function hart6(x)    
    alpha = [1.0 1.2 3.0 3.2]
    A = [10 3 17 3.5 1.7 8;
         0.05 10 17 0.1 8 14;
         3 3.5 1.7 10 17 8;
         17 8 0.05 10 0.1 14];
    P = 10^(-4) * [1312 1696 5569 124 8283 5886;
                   2329 4135 8307 3736 1004 9991;
                   2348 1451 3522 2883 3047 6650;
                   4047 8828 8732 5743 1091 381]
    
    outer = 0
    for ii = 1:4
        inner = 0
        for jj = 1:6
            xj = x[jj]
            Aij = A[ii, jj]
            Pij = P[ii, jj]
            inner = inner + Aij*(xj-Pij)^2
        end
        new = alpha[ii] * exp(-inner)
        outer = outer + new
    end
    
    y = -outer
end



using LatinHypercubeSampling
using ScatteredInterpolation
using BlackBoxOptim
using Parameters
using Distances
using Statistics
using StatsBase
using LinearAlgebra
using Distributed
using StaticArrays
using NearestNeighbors

#Extended methods
import Base.minimum
import Base.maximum
import Statistics.std
import Statistics.median
import Statistics.mean

include("types.jl")
include("interface.jl")
include("rbf_hypers_opt.jl")
include("rbf_hypers_opt_utilities.jl")
include("sample_infill.jl")
include("scaling.jl")
include("smoptimize.jl")
include("search_range.jl")


#TODO
#####Tracing - Handle in nice way with option :verbose or :silent
#####Remove outliers from std infill. Gather statistics on this if it should be used
#####PlotlyJS parallel coordinates
#####Ability to retain sm params. Do this with an option or result type
#####Options for infill
#####Gather statistics, rippa or not? Then multiple regression points or not?
#####Statistics, infill criteria, choose different. Recursive criteria, like every other is distance and every other is min?
#####Statistics, infill criteria of min-std and min. Always take min and find min-std with the greatest distance from min.
#####keep adding N number of such points.
#####Result type containing
        # Original samples
        # Infill points
        # Smallest evaluated function value per iteration
        # Options?
        # Num iterations
        # ?

#####Infill
    #Verify that the infill works as expected with small test function
    #Handle several infill methods
    #Pareto front type of infill? 
    #How to define the distance between samples of a pareto front?
    #New designs.
    #1. Pure Exploration of design space. 
#Done   #1. Infill distance, find largest empty space - GA to find this, KD-tree to increase speed
        #2. Think of dimension scaling. Not possible with several interpolants?
#Done   #3. Find maximum STD
        #4. Large gradient
    #2. Pure Exploitation
#Done   #1. Take overall min prediction
#Done   #2. Take min median prediction
    #3. Both
#Done   #1. min(Median prediction - 1STD)
#Done   #2. min(Median prediction - 2STD)
        #3. min(Median prediction) >= c Distance from any other point
        #4. min(Median prediction - xSTD) where X is the smallest value achieving distance
        #   >= c from any other point

            
        

end # module
