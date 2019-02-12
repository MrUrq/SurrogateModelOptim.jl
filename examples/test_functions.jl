using Parameters

"""
    TestFunction

Data structure to store test functions.
"""
@with_kw struct TestFunction
    fun
    sr::Array{Tuple{Float64,Float64},1}
    min_val::Float64
    min_loc::Array{Float64,2}
    max_val::Float64
    max_loc::Array{Float64,2}
end

function rosenbrock_2D(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function rotatedHyperElipsoid_2D(x)
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

function styblinskiTang_2D(x)
    out = 0.0    
    xy = [x[1],x[2]]
    for i = 1:2
        out += xy[i]^4 - 16*xy[i]^2 + 5*xy[i]
    end
    out *= 0.5
    return out
end

function threeHumpCamel_2D(x)
    out = 2*x[1]^2 - 1.05*x[1]^4 + (x[1]^6)/6 + x[1]*x[2] + x[2]^2
end

function bohachevsky_2D(x)
    out = x[1]^2 + 2*x[2]^2 - 0.3*cos(3*π*x[1]) - 0.4*cos(4*π*x[2]) + 0.7
end

function permfunc_2D(x)
    out = 0.0
    xy = x
    for i = 1:2
        inner = 0.0
        for j = 1:2
            xj = xy[j]
            inner = inner + (j)*(xj^i-(1/j)^i)
        end
        out = out + inner^2
    end
    return out
end

function sphereFunction_2D(x)
    out = 0.0
    xy = x
    for i = 1:2
        out += xy[i]^2
    end
    return out
end

function sumOfDifferentPowers_2D(x)
    out = 0.0
    xy = x
    for i = 1:2
        out += abs(xy[i])^(i+1)
    end
    return out
end

function sumOfSquares_2D(x)
    out = 0.0
    xy = x 
    for i = 1:2
        out += i*xy[i]^2
    end
    return out
end

function trid_2D(x)
    xy = x
    out = (xy[1]-1)^2
    inner = 0.0    
    for i = 2:2
        out += (xy[i]-1)^2
        inner += xy[i]*xy[i-1]       
    end
    return out - inner
end

function booth_2D(x)
    out = (x[1]+2*x[2]-1)^2 + (2*x[1]+x[2]-5)^2 
end

function matyas_2D(x)
    out = 0.26*(x[1]^2+x[2]^2) - 0.48*x[1]*x[2]
end

function mccormick_2D(x)
    out = sin(x[1]+x[2]) + (x[1]-x[2])^2 - 1.5*x[1] + 2.5*x[2] + 1
end

function powersum_2D(x)
    out = 0.0
    b = [8.0,18.0]
    xy = x
    for i = 1:2
        inner = 0.0
        for j = 1:2
            inner += xy[j]^i
        end
        out += (inner-b[i])^2
    end
    return out
end

function zakharov_2D(x)
    out = 0.0
    tmp = 0.0
    xy = x
    for i = 1:2
        out += xy[i]^2
        tmp += 0.5*i*xy[i]
    end
    return out + tmp^2 + tmp^4
end

function sixHumpCamel_2D(x)
    out = (4-2.1*x[1]^2 + (x[1]^4)/3)*x[1]^2 + x[1]*x[2] + (-4 + 4*x[2]^2)*x[2]^2
end

function dixonPrice_2D(x)
    out = 0.0
    xy = x
    out += (x[1]-1)^2
    for i = 2:2
        out += i*(2*xy[i]^2-xy[i]-1)^2
    end
    return out
end

function beale_2D(x)
    out = (1.5-x[1]+x[1]*x[2])^2 + (2.25-x[1]+x[1]*x[2]^2)^2 + (2.625-x[1]+x[1]*x[2]^3)^2
end

function branin_2D(x)
    out = 1*(x[2] - 5.1/(4*π^2)*x[1]^2 + 5/π*x[1] - 6)^2 + 10*(1-1/8/π)*cos(x[1]) + 10
end

function goldsteinPrice_2D(x)
    out = (1 + ((x[1]+x[2]+1)^2)*(19-14*x[1]+3*x[1]^2-14*x[2]+6*x[1]*x[2]+3*x[2]^2)) * 
                      (30 + ((2*x[1]-3*x[2])^2)*(18-32*x[1]+12*x[1]^2+48*x[2]-36*x[1]*x[2]+27x[2]^2))
end

function permdbeta_2D(x)
    out = 0.0
    xy = x    
    for i = 1:2
        inner = 0.0
        for j = 1:2
            inner += (j^i+0.5)*((xy[j]/j)^i-1)
        end
        out += inner^2
    end
    return out    
end

function rastrigin_2D(x)
    sum = 0
    for ii = 1:2
        xi = x[ii]
        sum = sum + (xi^2 - 10*cos(2*pi*xi))
    end
    out = 10*2 + sum
end

function hart_4D(x)

    alpha = [1.0 1.2 3.0 3.2]
    A = [10 3 17 3.5 1.7 8;
        0.05 10 17 0.1 8 14;
        3 3.5 1.7 10 17 8;
        17 8 0.05 10 0.1 14]
    P = 10^(-4) * [1312 1696 5569 124 8283 5886;
                2329 4135 8307 3736 1004 9991;
                2348 1451 3522 2883 3047 6650;
                4047 8828 8732 5743 1091 381]

    outer = 0
    for ii = 1:4
        inner = 0
        for jj = 1:4
            xj = x[jj]
            Aij = A[ii, jj]
            Pij = P[ii, jj]
            inner += Aij*(xj-Pij)^2
        end
        outer += alpha[ii] * exp(-inner)
    end

    out = (1.1 - outer) / 0.839

end

function hart_6D(x)    
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

function rosenbrock_ND(x)
    d = length(x)
    sum = 0.0
    for ii = 1:(d-1)
        xi = x[ii]
        xnext = x[ii+1]
        sum += 100*(xnext-xi^2)^2 + (xi-1)^2        
    end
    out = sum
end

##### Note, the maximums are approximate. The functions may have several optima, only one is stored.
test_funs = Dict(   
    :rosenbrock_2D => TestFunction( fun = rosenbrock_2D,
                                    sr = [(-5.0,5.0),(-5.0,5.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([1.0; 1.0]'),
                                    max_val = 90036.0,
                                    max_loc = permutedims([-5.0; -5.0]')),

    :rotatedHyperElipsoid_2D => TestFunction( fun = rotatedHyperElipsoid_2D,
                                    sr = [(-65.536, 65.536), (-65.536, 65.536)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 12884.901888,
                                    max_loc = permutedims([65.536; -65.536]')),

    :styblinskiTang_2D => TestFunction( fun = styblinskiTang_2D,
                                    sr = [(-5.0, 5.0), (-5.0, 5.0)],
                                    min_val = -78.33233140754285,
                                    min_loc = permutedims([-2.90353; -2.90353]'),
                                    max_val = 250.0,
                                    max_loc = permutedims([5.0; 5.0]')),    

    :threeHumpCamel_2D => TestFunction( fun = threeHumpCamel_2D,
                                    sr = [(-5.0, 5.0), (-5.0, 5.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 2047.9166666666665,
                                    max_loc = permutedims([5.0; 5.0]')),       
                                    
    :bohachevsky_2D => TestFunction( fun = bohachevsky_2D,
                                    sr = [(-100.0, 100.0), (-100.0, 100.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 30000.0,
                                    max_loc = permutedims([-100.0; -100.0]')),   
    
    :permfunc_2D => TestFunction( fun = permfunc_2D,
                                    sr = [(-2.0, 2.0), (-2.0, 2.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([1.0; 0.5]'),
                                    max_val = 174.25,
                                    max_loc = permutedims([-2.0; -2.0]')),

    :sphereFunction_2D => TestFunction( fun = sphereFunction_2D,
                                    sr = [(-5.12, 5.12), (-5.12, 5.12)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 52.4288,
                                    max_loc = permutedims([5.12; 5.12]')),

    :sumOfDifferentPowers_2D => TestFunction( fun = sumOfDifferentPowers_2D,
                                    sr = [(-1.0, 1.0), (-1.0, 1.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 2.0,
                                    max_loc = permutedims([1.0; -1.0]')),
    
    :sumOfSquares_2D => TestFunction( fun = sumOfSquares_2D,
                                    sr = [(-10.0, 10.0), (-10.0, 10.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0;0.0]'),
                                    max_val = 300.0,
                                    max_loc = permutedims([-10.0; 10.0]')),

    :trid_2D => TestFunction( fun = trid_2D,
                                    sr = [(-4.0, 4.0), (-4.0, 4.0)],
                                    min_val = -2.0,
                                    min_loc = permutedims([2.0; 2.0]'),
                                    max_val = 50.0,
                                    max_loc = permutedims([-4.0; 4.0]')),

    :booth_2D => TestFunction( fun = booth_2D,
                                    sr = [(-10.0, 10.0), (-10.0, 10.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([3.0; -1.0]'),
                                    max_val = 2186.0,
                                    max_loc = permutedims([-10.0; -10.0]')),

    :matyas_2D => TestFunction( fun = matyas_2D,
                                    sr = [(-10.0, 10.0), (-10.0, 10.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 100.0,
                                    max_loc = permutedims([10.0; -10.0]')),

    :mccormick_2D => TestFunction( fun = mccormick_2D,
                                    sr = [(-1.5, 4.0), (-3.0, 4.0)],
                                    min_val = -1.9132229549810367,
                                    min_loc = permutedims([-0.547198; -1.5472]'),
                                    max_val = 44.09847214410396,
                                    max_loc = permutedims([-1.5; 4.0]')),

    :powersum_2D => TestFunction( fun = powersum_2D,
                                    sr = [(0.0, 2.0), (-0.0, 2.0)],
                                    min_val = 116.0,
                                    min_loc = permutedims([2.0; 2.0]'),
                                    max_val = 388.0,
                                    max_loc = permutedims([0.0;0.0]')),

    :zakharov_2D => TestFunction( fun = zakharov_2D,
                                    sr = [(-5.0, 10.0), (-5.0, 10.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([0.0; 0.0]'),
                                    max_val = 51050.0,
                                    max_loc = permutedims([10.0; 10.0]')),

    :sixHumpCamel_2D => TestFunction( fun = sixHumpCamel_2D,
                                    sr = [(-3.0, 3.0), (-2.0, 2.0)],
                                    min_val = -1.0316284534898774,
                                    min_loc = permutedims([-0.089842; 0.712656]'),
                                    max_val = 162.89999999999998,
                                    max_loc = permutedims([3.0; 2.0]')),

    :dixonPrice_2D => TestFunction( fun = dixonPrice_2D,
                                    sr = [(-10.0, 10.0), (-10.0, 10.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([1.0; -0.5]'),
                                    max_val = 87483.0,
                                    max_loc = permutedims([-10.0; -10.0]')),

    :beale_2D => TestFunction( fun = beale_2D,
                                    sr = [(-4.5, 4.5), (-4.5, 4.5)],
                                    min_val = 0.0,
                                    min_loc = permutedims([3.0; -0.5]'),
                                    max_val = 181853.61328125,
                                    max_loc = permutedims([-4.5; -4.5]')),

    :branin_2D => TestFunction( fun = branin_2D,
                                    sr = [(-5.0, 10.0), (0.0, 15.0)],
                                    min_val = 0.39788735772973816,
                                    min_loc = permutedims([3.14159; 2.275]'),
                                    max_val = 308.12909601160663,
                                    max_loc = permutedims([-5.0; 1.12981e-15]')),

    :goldsteinPrice_2D => TestFunction( fun = goldsteinPrice_2D,
                                    sr = [(-2.0, 2.0), (-2.0, 2.0)],
                                    min_val = 3,
                                    min_loc = permutedims([0; -1.0]'),
                                    max_val = 1.0156902717980599e6,
                                    max_loc = permutedims([-1.73737; 2.0]')),

    :permdbeta_2D => TestFunction( fun = permdbeta_2D,
                                    sr = [(-2.0, 2.0), (-2.0, 2.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([1.0; 2.0]'),
                                    max_val = 110.5,
                                    max_loc = permutedims([-2.0; -2.0]')),

    :rastrigin_2D => TestFunction( fun = rastrigin_2D,
                                    sr = [(-5.12, 5.12), (-5.12, 5.12)],
                                    min_val = 0.0,
                                    min_loc = permutedims([1.05917e-9; 1.72311e-9]'),
                                    max_val = 80.70658038767792,
                                    max_loc = permutedims([-4.52299; 4.52299]')),

    :hart_4D => TestFunction( fun = hart_4D,
                                    sr = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
                                    min_val = -3.1344941412224,
                                    min_loc = permutedims([0.187395; 0.194152; 0.557918; 0.26478]'),
                                    max_val = 1.3095406214559708,
                                    max_loc = permutedims([1.0; 1.0; 1.42362e-13; 1.0]')),

    :hart_6D => TestFunction( fun = hart_6D,
                                    sr = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
                                    min_val = -3.322368011415515,
                                    min_loc = permutedims([0.20169; 0.150011; 0.476874; 0.275332; 0.311652; 0.657301]'),
                                    max_val = -2.8124505439686604e-8,
                                    max_loc = permutedims([1.0; 1.0; 3.13887e-17; 1.0; 1.0; 1.0]')),

    :rosenbrock_12D => TestFunction( fun = rosenbrock_ND,
                                    sr = [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)],
                                    min_val = 0.0,
                                    min_loc = permutedims([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]'),
                                    max_val = 990396.0,
                                    max_loc = permutedims([-5.0; -5.0; -5.0; -5.0; -5.0; -5.0; -5.0; -5.0; -5.0; -5.0; -5.0; -5.0]')),
    
                                    );