
# # D dispersion op
# # N nonlinear term
# using FFTW

# function integrate(u₀, L)

# 	dz = L / 1000
# 	z = 0
# 	l = []
# 	k₅ = nothing

# 	D = 1.0 + zeros(length(u₀))

# 	N(û) = fft(ifft())

# 	while z < L
# 		û, err = _erk4ip_step!(û)
# 		z = z + dz

# 		append!(l, z)
# 	end


# 	function _erk4ip_step!(û, k₅ = nothing)


# 		a_int = dispersion_op * û
# 		if isnothing(k₅)
# 			k₁ = dz * D * N(û)
# 		else
# 			k₁ = D * k₅
# 		end
# 		k₂ = dz * N(a_int + k₁ / 2)
# 		k₃ = dz * N(a_int + k₃ / 2)
# 		k₄ = dz * N(D * (a_int + k₃))
# 		β = D * (a_int + k₁ / 6 + k₂ / 3 + k₃ / 3)
# 		sol₄ = β + k₄ / 6
# 		k₅ = dz * N(sol₄)
# 		sol₃ = beta + k₄ / 15 + k₅ / 10
# 		error = _compute_error(sol₄, sol₃)
# 		return sol₄, error


# 	end

# 	function _compute_error(a, b)
# 		sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
# 	end

# end


function _erk4ip_step!(stepper, model, z, reltol, maxiters)
	e = exp.(0.5 * stepper.dz * model.dispersion_term)

	Uip = e .* stepper.U

	k1 = e .* stepper.NU

	k2 = model.nonlinear_function(model.ifftp * (Uip .+ 0.5 * stepper.dz * k1), model)

	k3 = model.nonlinear_function(model.ifftp * (Uip .+ 0.5 * stepper.dz * k2), model)

	k4 = model.nonlinear_function(model.ifftp * (e .* (Uip .+ stepper.dz * k3)), model)

	r = e .* (Uip .+ stepper.dz * (k1 / 6.0 .+ k2 / 3.0 .+ k3 / 3.0))

	U1 = r .+ stepper.dz * k4 / 6.0

	stepper.k5 = model.nonlinear_function(model.ifftp * U1, model)

	U2 = r .+ stepper.dz * (k4 / 15.0 .+ stepper.k5 / 10.0)

	stepper.local_error = _compute_error(U1, U2)

	dzopt =
		max(0.5, min(2.0, 0.9 * sqrt(sqrt(reltol / stepper.local_error)))) * stepper.dz

	if stepper.local_error <= reltol

		stepper.dz = min(dzopt, abs(z - stepper.z))
		stepper.z = stepper.z + stepper.dz
		stepper.U = U1
		stepper.NU = stepper.k5
	else
		stepper.dz = dzopt
		stepper.it = stepper.it + 1
		if (stepper.it >= maxiters)
			throw(ErrorException("Max number of iteration exceeded!"))
		end
	end
end
