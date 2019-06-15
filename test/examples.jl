@testset "examples" begin
    @test include("../examples/test_function_optimisation.jl")
    @test include("../examples/test_function_optimisation_noise.jl")
    @test include("../examples/test_function_optimisation_categorical.jl")
    @test include("../examples/test_function_optimisation_discontinuous.jl")
end