# PROJECT_POS    Project a point onto the nonnegative orthant.
#
#   project_pos(v) is the positive part of v.
function project_pos(v)
    x = max(v, 0)
    return x
end
