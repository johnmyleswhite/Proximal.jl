# PROX_SEPARABLE   Evaluate the prox operator of a fully separable function.
#
# Arguments:
#
#  v is the point at which to evaluate the operator.
#  fp is a subgradient oracle for the function.
#  lambda (optional) is the proximal parameter; defaults to 1.
#  l (optional) is a lower bound for x; defaults to -Inf.
#  u (optional) is an upper bound for x; defaults to Inf.
#  x0 (optional) is a value at which to warm start the algorithm.
#  tol (optional) is a stopping tolerance.
#  MAX_ITER (optional) is the maximum number of iterations.
#
# Examples:
#
#  v = randn(n,1);
#  x = prox_separable(v, 1, @(w) sign(w));
#  [x iter] = prox_separable(v, 1, @(w) sign(w));
#
# This function can be called in vectorized form if fp is vectorized,
# i.e., if fp works elementwise, then v, l, u, and x0 can be vectors.
function prox_separable(
    v,
    lambda::Real = NaN,
    fp::Function,
    l,
    u,
    x0,
    tol::Real = NaN,
    maxiter::Integer = 0
)
    n = length(v)

    if isnan(lambda)
        lambda = 1.0
    end
    rho = 1 / lambda

    if any(isnan(l)) || isempty(l)
        l = -inf(n, 1)
    end
    if any(isnan(u)) || isempty(u)
        u = inf(n, 1)
    end
    if any(isnan(x0)) || isempty(x0)
        x0 = zeros(n, 1)
    end
    if isnan(tol)
        tol = 1e-8
    end
    if maxiter == 0
        maxiter = 500
    end

    iter = 0
    x = max(l, min(x0, u))

    while any(u .- l .> tol) && iter < MAX_ITER
        g = fp(x) .+ rho .* (x .- v)

        idx = g .> 0
        l[idx] = max(l[idx], x[idx] .- g[idx] ./ rho)
        u[idx] = x[idx]

        idx = !idx
        u[idx] = min(u[idx], x[idx] .- g[idx] ./ rho)
        l[idx] = x[idx]

        x = (l .+ u) ./ 2
        iter += 1
    end

    if any(u .- l .> tol)
        @printf(
            STDERR,
            "Warning: %d entries did not converge; max interval size = %f.",
            sum(u .- l .> tol),
            max(u .- l)
        )
    end

    return x, iter
end
