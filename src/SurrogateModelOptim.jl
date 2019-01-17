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
using Printf

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

#####Fix minimum infill, should be median

#####Result type containing
        # Original samples
        # Infill points
        # Smallest evaluated function value per iteration
        # Options?
        # Num iterations
        # Interpolant hypers
#####Fancy results printing

#####Go through test functions find max, min, range §and create a type that holds all the info

#####Gather statistics, rippa or not? 
#####Multiple regression points or not?

#####PlotlyJS parallel coordinates
#####Plot accuracy as a function of "time" predicted and actual

#####Performance of several infill points at once vs one at the time?

#####Outlier removal z from std infill. Gather statistics on this if it should be used
#####Statistics, infill criteria, choose different. Recursive criteria, like every other is distance and every other is min?
#####Statistics, infill criteria of min-std and min. Always take min and find min-std with the greatest distance from min.

#####Infill
    #Verify that the infill works as expected with small test function. 
    #Plot it and show the where the points are added after each iteration
    #Pareto front infill option possibly?
    
    #1. Pure Exploration of design space. 
#Done   #1. Infill distance, find largest empty space - GA to find this, KD-tree to increase speed
        #2. Think of dimension scaling. Not possible with several interpolants?   - Outside of project scope
#Done   #3. Find maximum STD
    #2. Pure Exploitation
#Done   #1. Take overall min prediction
#Done   #2. Take min median prediction
    #3. Both
#Done   #1. min(Median prediction - 1STD)
#Done   #2. min(Median prediction - 2STD)
        #3. min(Median prediction) >= c Distance from any other point
        #4. min(Median prediction - xSTD) where X is the smallest value achieving distance
        #   >= c from any other point

#####Integer test case
        

end # module
