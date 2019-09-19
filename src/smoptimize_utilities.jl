"""
    f_opt_eval(f,plan,samples;trace::Symbol=:silent)

Calculate the objective function value(s) and provide intermediate results printing
showing improvements over previous best iteration.
"""
function f_opt_eval(f,plan,samples;trace::Symbol=:silent)

    (trace == :verbose) && println("Evaluating function ",size(plan,2)," times")

    new_samples = mapslices(f,plan,dims=1) 

    print_f_opt(trace,new_samples,samples,plan)

    return new_samples
end

function print_f_opt(trace,new_samples,samples,plan)

    new_min = minimum(new_samples)
    old_min = minimum(samples)
    new_max = maximum(new_samples)
    old_max = maximum(samples)

    if trace == :compact
        min_sample = minimum((old_min,new_min))
        if isless(new_min,old_min)
            print(@sprintf("%20.7g %20.7g\n", min_sample, old_min-new_min))
        else
            print(@sprintf("%20.7g %20.7s\n", min_sample, "N/A"))
        end
        
    elseif trace == :verbose
        println("    ---------------------------------------------------------------")
        print("    ")
        for sample in new_samples
            printstyled(@sprintf("%-15.7g",sample))
        end
        println("\t function evaluation\n")

        print("    Minimum sample value ")
        printstyled(@sprintf("%.7g",min(new_min,old_min)); color=:light_green, bold=true)
        print("\t(maximum = ")
        printstyled(@sprintf("%.7g",max(new_max,old_max));)
        print(")")

        if isless(new_min,old_min)
            print("\t (Improvement from last iteration ", @sprintf("%.7g",old_min-new_min),")\n")
        else
            print("\t (Improvement from last iteration N/A)\n")
        end
        print("\n")
    end
end

"""
    f_opt_eval(f,plan;trace::Symbol=:silent)

Calculate the objective function value(s) and plot value if trace.
"""
function f_opt_eval(f,plan;trace::Symbol=:silent)

    (trace == :verbose) && println("Evaluating function ",size(plan,2)," times")
    (trace == :compact) && print(@sprintf("%s %21s %20s\n", "   Iteration", "Function value", "Improvement"))

    new_samples = mapslices(f,plan,dims=1) 

    if trace == :verbose
        new_min = minimum(new_samples)
        new_max = maximum(new_samples)

        print("    Minimum sample value ")
        printstyled(@sprintf("%.7g",new_min); color=:light_green, bold=true)
        print("\t(maximum = ")
        printstyled(@sprintf("%.7g",new_max);)
        println(")\n")
    end

    return new_samples
end

"""
    print_iteration(trace,i,iterations)

REPL printing of the current iteration if trace.
"""
function print_iteration(trace,i,iterations)
    if trace == :compact
        print(@sprintf("%4.4g%s%g",i, " out of ", iterations))
    end

    if trace == :verbose
        print("\n\nFunction eval. ")
        printstyled("iteration ";bold=true,color=:light_red)
        printstyled(i,bold=true)
        print(" out of ", iterations)
        print("\n")
    end
end

"""
    best_fitness(result::SurrogateResult)

Extract the best fitness from the `SurrogateResult` type.
"""
function best_fitness(result::SurrogateResult)
    minimum([result.lhc_samples result.infill_samples])
end

"""
    best_candidate(result::SurrogateResult)

Extract the design parameters yielding the best fitness from the `SurrogateResult` type.
"""
function best_candidate(result::SurrogateResult)
    min_ind = findmin([result.lhc_samples result.infill_samples])[2][2]
    return [result.lhc_plan result.infill_plan][:,min_ind:min_ind]
end

"""
    f_calls(result::SurrogateResult)

Extract the number of function calls used from the `SurrogateResult` type.
"""
function f_calls(result::SurrogateResult)
    return length([result.lhc_samples result.infill_samples])
end