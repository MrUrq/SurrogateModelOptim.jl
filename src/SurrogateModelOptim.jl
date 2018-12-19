module SurrogateModelOptim

export  smoptimize,
        rosenbrock2d

function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

using LatinHypercubeSampling: LHCoptim
using ScatteredInterpolation
using BlackBoxOptim: bboptimize
using Parameters: reconstruct, type2dict, @with_kw, @unpack
using Distances
using Statistics
using LinearAlgebra
using Distributed
import Base.minimum


include("types.jl")
include("interface.jl")
include("rbf_hypers_opt.jl")
include("sample_infill.jl")
include("scaling.jl")
include("smoptimize.jl")


end # module
