# PROX_MAX    The proximal operator of the max function.
#
#   prox_max(v,lambda) is the proximal operator of the max
#   of the entries of v. This function is not vectorized.

function prox_max(v, lambda)
    tol = 1e-8
    maxiter = 100
    rho = 1 / lambda

    n = length(v)
    tl = min(v) - 1 / n
    tu = max(v)

    function g(t)
        return sum(max(0, rho * (v - t))) - 1
    end

    iter = 0
    while tu .- tl .> tol && iter < maxiter
        t0 = (tl + tu) / 2
        if sign(g[t0]) == sign(g[tl])
            tl = t0
        else
            tu = t0
        end 
        iter += 1
    end

    x = min(t0, v)

    if tu .- tl .> tol
        warn("Algorithm did not converge.")
    end

    return x
end
