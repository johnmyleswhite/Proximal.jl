# PROJECT_AFFINE    Project a point into an affine set.
#
#   project_affine(v,A,b) is the projection of v onto
#   the affine set { x | Ax = b }.
#
#   You can also call the function as
#
#     [x C d] = project_affine(v,A,b);
#
#   and then call it again with different argument v2 via
#
#     x2 = project_affine(v2,A,b,C,d);
#
#   If calling the function repeatedly with the same A and b,
#   all evaluations after the initial one will be much faster
#   when passing in the cached values C and d.

# TODO: Determine dimensions of all objects
function project_affine(
	v::Array,
	A::Array,
	b::Array,
	C::Array = Array(Float64, 0),
	d::Array = Array(Float64, 0)
)
	# x = v - pinv(A) * (A * v - b)
    if isempty(C) || isempty(d)
        pA = pinv(A)
        C = eye(length(v)) - pA * A
        d = pA * b
    end

	x = C * v + d

	return x, C, d
end
