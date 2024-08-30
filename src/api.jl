




function GNLSEProblem(t, wg::Waveguide)

	fftp = plan_fft(t, flags = FFTW.MEASURE)
	ifftp = plan_ifft(t, flags = FFTW.MEASURE)

	T = t[end] - t[1]
	N = length(t)
	dt = T / N

	ν = fftshift((-N÷2:N÷2-1) ./ T)
	ω = 2π .* ν

	ω0 = 2pi * c / wg.λc

	# default values if no raman 
	raman_freq_response = nothing
	# raman prob data
	if wg.raman_model.fr != 0.0
		raman_freq_response = conj((fftp * ifftshift(wg.raman_model.time_response(t))))
	end

	nonlinear_function = choose_nonlinear_term(wg.self_steepening, !isnothing(raman_freq_response))

	dispersion_term = -0.5wg.α .+ 1im * sum([(wg.βs[i] / factorial(i)) .* ω .^ i for i in eachindex(wg.βs)])
	GNLSEProblem(ω, dt, N, fftp, ifftp, dispersion_term, nonlinear_function, wg.raman_model.fr, wg.γ, raman_freq_response, ω0, wg.L)
end




function _compute_error(a, b)
	sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
end

function _integrate_to_z(stepper, z, prob, maxiters, reltol)
	stepper.it = 0
	while stepper.z < z
		_erk4ip_step!(stepper, prob)
		stepper.local_error = _compute_error(stepper.U1, stepper.U2)

		dzopt =
			max(0.5, min(2.0, 0.9 * sqrt(sqrt(reltol / stepper.local_error)))) * stepper.dz

		if stepper.local_error <= reltol
			stepper.dz = min(dzopt, abs(z - stepper.z))
			stepper.z = stepper.z + stepper.dz
			stepper.U = stepper.U1
			stepper.NU = stepper.k5
		else
			stepper.dz = dzopt
			stepper.it = stepper.it + 1
			if (stepper.it >= maxiters)
				throw(ErrorException("Max number of iteration exceeded!"))
			end
		end
	end
end

function gnlse(u, t, prob::GNLSEProblem; nsaves = 20, dz = 1.0, reltol = 1e-6, maxiters = 1000)

	dz = min(prob.L / (2 * nsaves), dz)
	# initial stepper
	stepper = Stepper(prob.fftp * ComplexF64.(u), prob.nonlinear_function(u, prob), dz, 0.0, 0.0, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, 0)
	zsaves = (0:nsaves) * prob.L / nsaves
	M = zeros(ComplexF64, (nsaves + 1, prob.N))
	M[1, :] = stepper.U
	ϵ_hist = zeros(nsaves + 1)
	ϵ_hist[1] = 0.0

	for i ∈ 2:nsaves+1
		_integrate_to_z(stepper, zsaves[i], prob, maxiters, reltol)
		ϵ_hist[i] = stepper.local_error
		M[i, :] = stepper.U
	end

	return Solution(zsaves, t, prob.ω / 2pi, ifft(M, 2), M)
end

function gnlse(u, t, wg::Waveguide; args...)
	gnlse(u, t, GNLSEProblem(t, wg); args...)
end

function gnlse(u, t, probs::Vector{GNLSEProblem}; args...)
	sols = [gnlse(u, t, probs[1]; args...)]
	for prob ∈ probs[2:end]
		push!(sols, gnlse(sols[end].At[end, :], t, prob; args...))
	end
	combine(sols)
end

function gnlse(u, t, wgs::Vector{Waveguide}; args...)
	probs = [GNLSEProblem(t, wg) for wg in wgs]
	gnlse(u, t, probs; args...)
end
