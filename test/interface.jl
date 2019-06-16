@testset "interface" begin
    @test typeof(SurrogateModelOptim.Options()) == SurrogateModelOptim.Options
    @test_throws MethodError SurrogateModelOptim.Options(not_a_valid_arugment=true)
    @test_throws AssertionError SurrogateModelOptim.Options(smooth=true)
    @test typeof(SurrogateModelOptim.Options(smooth=false)) == SurrogateModelOptim.Options
    @test typeof(SurrogateModelOptim.Options(smooth=:variable)) == SurrogateModelOptim.Options
    @test typeof(SurrogateModelOptim.Options(smooth=:single)) == SurrogateModelOptim.Options
    @test typeof(SurrogateModelOptim.Options(smooth=:single_user)) == SurrogateModelOptim.Options
end