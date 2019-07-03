search_range = [(-10.0,10.0), (0.0,2.0), (25.0,40.0)]
n_samples = 6
plan = SurrogateModelOptim.scaled_LHC_sampling_plan(search_range,n_samples,1000;trace=:silent)

@testset "scaled LHC plan" begin
    @test vec(extrema(plan,dims=2)) == search_range
end

@testset "function eval" begin
    @test size(SurrogateModelOptim.f_opt_eval(minimum,plan)) == (1,n_samples)
    samples = mapslices(minimum,plan,dims=1)
    @test size(SurrogateModelOptim.f_opt_eval(minimum,plan,samples)) == (1,n_samples)
end