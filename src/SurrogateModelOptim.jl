module SurrogateModelOptim

export  smoptimize,
        model_infill,
        TestFunction


using LatinHypercubeSampling
using ScatteredInterpolation
using BlackBoxOptim
using Distributions
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


#####Result type containing
        # Original samples
        # Infill points
        # Smallest evaluated function value per iteration
        # Options?
        # Num iterations
        # Interpolant hypers
#####Fancy results printing
#####Call the method with the results type already?

#####Change distance to negative instead of 1/x? same with std?

#####Try the updated BBoptim and answer github thingy

#####Call the method with a test function + options shortcut

#####Create plot of optim method vs optim method, shaded area thing. Simplex?

#####Performance as a function of iteration number?

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
#####Documentation
#####Examples
#####Function headers

end # module
