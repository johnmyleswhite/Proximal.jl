# PROX_PRECOMPOSE    Proximal precomposition rule.
#
#   Let g(x) = f(tx + b). Then
#
#     prox_g = prox_precompose(prox_f,t,b).
#
#   Here, prox_f only takes v and lambda as arguments.
function prox_precompose(prox_f, t, b)
    function prox_g(v, lambda)
    	((1 / t) * (prox_f(v * t + b, 1 / (lambda * t^2)) - b))
    end
    return prox_g
end
