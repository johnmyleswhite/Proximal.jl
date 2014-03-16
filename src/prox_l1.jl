# PROX_L1    The proximal operator of the l1 norm.
#
#   prox_l1(v,lambda) is the proximal operator of the l1 norm
#   with parameter lambda.
function x = prox_l1(v, lambda)
    x = max(0, v - lambda) - max(0, -v - lambda)
    return x
end
