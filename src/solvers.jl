

# function _erk4ip_step!(û, k₅ = nothing)

# 	a_int = dispersion_op * û
# 	if isnothing(k₅)
# 		k₁ = dz * D * N(û)
# 	else
# 		k₁ = D * k₅
# 	end
# 	k₂ = dz * N(a_int + k₁ / 2)
# 	k₃ = dz * N(a_int + k₃ / 2)
# 	k₄ = dz * N(D * (a_int + k₃))
# 	β = D * (a_int + k₁ / 6 + k₂ / 3 + k₃ / 3)
# 	sol₄ = β + k₄ / 6
# 	k₅ = dz * N(sol₄)
# 	sol₃ = beta + k₄ / 15 + k₅ / 10
# 	error = _compute_error(sol₄, sol₃)
# 	return sol₄, error

# end
