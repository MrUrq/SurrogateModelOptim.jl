@testset "surrogate" begin
   
    function rosenbrock_2D(x)
        return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    end

    # Try whole run
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=1,
                        num_start_samples=2,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = true,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test typeof(result) == SurrogateModelOptim.SurrogateResult

    # Try with more than 1 interpolants
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=2,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = true,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test length(result.sm_interpolant_settings) == 2

    # Try with more than 2 start samples
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = true,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test length(result.lhc_samples) == 5

    # Try with more than 1 infill points
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=2,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = true,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test length(result.infill_samples) == 3*2 #iterations*num_interpolants

    # Try that it works without using rippas algorithm
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=2,
                        trace=false,
                        rippa = false,
                        variable_kernel_width = true,
                        variable_dim_scaling = true,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test typeof(result) == SurrogateModelOptim.SurrogateResult

    # Try that it works without variable kernel hyperparameter optimisation
    #First try that it works with variable kernel optimisation
    @test typeof(result.sm_interpolant_settings[1].kernelFunc) == Array{SurrogateModelOptim.ScatteredInterpolation.RadialBasisFunction,1}
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=2,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = false,
                        variable_dim_scaling = true,
                        kerns = [SurrogateModelOptim.ScatteredInterpolation.Gaussian],
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test typeof(result.sm_interpolant_settings[1].kernelFunc) == SurrogateModelOptim.ScatteredInterpolation.Gaussian{Float64}


    # Try without variable dim scaling
    # First try that it works as intended with dim scaling
    @test typeof(result.sm_interpolant_settings[1].scaling) == Array{Float64,1}
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = false,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test result.sm_interpolant_settings[1].scaling == false

    # Try smoothing off
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = false,
                        smooth=:false,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test result.sm_interpolant_settings[1].smooth == false

    # Try smoothing variable for each point
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = false,
                        smooth=:variable,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test typeof(result.sm_interpolant_settings[1].smooth) == Array{Float64,1}

    # Try smoothing with a single value for all points
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = false,
                        smooth=:single,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test typeof(result.sm_interpolant_settings[1].smooth) == Float64

    # Try settings the smoothing from the interface
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=3,
                        num_interpolants=2,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = false,
                        smooth=:single_user,
                        smooth_user = 1337.0,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test result.sm_interpolant_settings[1].smooth == 1337.0

    # Test that the infill functions work
    result = smoptimize(rosenbrock_2D, [(-5.0,5.0),(-5.0,5.0)];
                        options=SurrogateModelOptim.Options(
                        iterations=6,
                        num_interpolants=1,
                        num_start_samples=5,
                        rbf_opt_gens=50,
                        infill_iterations=50,
                        num_infill_points=1,
                        trace=false,
                        rippa = true,
                        variable_kernel_width = true,
                        variable_dim_scaling = false,
                        parallel_surrogate=false,
                        infill_funcs=[:min,:mean,:median,:std,:dist]
                            ));
    @test result.infill_type == [:min,:mean,:median,:std,:dist,:min]
end