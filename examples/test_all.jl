dir = @__DIR__
include(joinpath(dir,"test_function_optimisation.jl"))
include(joinpath(dir,"test_function_optimisation_noise.jl"))
include(joinpath(dir,"test_function_optimisation_multi_obj.jl"))
include(joinpath(dir,"test_function_optimisation_fast_options.jl"))
include(joinpath(dir,"test_function_optimisation_discontinuous.jl"))
include(joinpath(dir,"test_function_optimisation_constraint.jl"))
include(joinpath(dir,"test_function_optimisation_categorical.jl"))