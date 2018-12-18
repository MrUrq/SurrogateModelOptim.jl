using SurrogateModelOptim
using Test

@testset "Interface" begin
    @testset "Default types" begin
        #Wrong use of interface
        @test_throws AssertionError SurrogateModelOptim._setup_parameters((-1,1.0))
        @test_throws AssertionError SurrogateModelOptim._setup_parameters((-1.0,1.0))
        @test_throws AssertionError SurrogateModelOptim._setup_parameters([(-1.0,1.0)], 1)

        #Correct use of interface
        @test typeof(SurrogateModelOptim._setup_parameters((-1.0,1.0),2)) == 
            SurrogateModelOptim._Parameters
        @test typeof(SurrogateModelOptim._setup_parameters([(-1.0,1.0)])) ==
            SurrogateModelOptim._Parameters
    end

    @testset "Rule of thumb initial sample number" begin
        n_samples = 3
        rule_of_thumb = 10
        @test SurrogateModelOptim._setup_parameters((-1.0,1.0),n_samples).NumStartSamples == n_samples*rule_of_thumb
        @test SurrogateModelOptim._setup_parameters(repeat([(-1.0,1.0)],n_samples)).NumStartSamples == n_samples*rule_of_thumb
    end
end

@testset "Sampling plan" begin
    SearchRange = (-2.0,1.0)
    NumDimensions = 2
    NumStartSamples = 5
    SamplingPlanOptGens = 100
    @test extrema(SurrogateModelOptim._LHC_sampling_plan(SearchRange,NumDimensions,NumStartSamples,SamplingPlanOptGens)) == SearchRange


    SearchRange = [(-2.0,1.0),(-3.0,0.5)]
    NumDimensions = false
    NumStartSamples = 5
    SamplingPlanOptGens = 100
    @test extrema(SurrogateModelOptim._LHC_sampling_plan(SearchRange,NumDimensions,NumStartSamples,SamplingPlanOptGens)[1,:]) == SearchRange[1]
    @test extrema(SurrogateModelOptim._LHC_sampling_plan(SearchRange,NumDimensions,NumStartSamples,SamplingPlanOptGens)[2,:]) == SearchRange[2]
end