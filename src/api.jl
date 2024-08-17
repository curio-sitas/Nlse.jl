

function choose_nonlinear_term(self_steepening::Bool = true, raman::Bool = true)

	if self_steepening && raman
		return NL_spm_self_steepening_raman
	elseif self_steepening && !raman
		return NL_spm_self_steepening
	elseif !self_steepening && raman
		return NL_spm_raman
	else
		return NL_spm
	end

end



function create_model(u::AbstractArray{ComplexF64}, t, wg::Waveguide; self_steepening::Bool = false, raman::Bool = false, normalize = true)

	fftp = plan_fft(u)
	ifftp = plan_ifft(u)

	T = t[end] - t[1]
	N = length(u)
	dt = T / N

	ν = fftshift((-N÷2:N÷2-1) ./ T)
	ω = 2π .* ν

	ω0 = 2pi * c / wg.λc

	# default values if no raman 
	raman_freq_response = nothing
	# raman model data
	if raman && wg.raman_model.fr != 0.0
		raman_freq_response = conj((fftp * ifftshift(wg.raman_model.time_response(t))))
	end

	dispersion_term = -0.5wg.α .+ 1im * sum([(wg.βs[i] / factorial(i)) .* ω .^ i for i in eachindex(wg.βs)])
	nonlinear_function = choose_nonlinear_term(self_steepening, raman)
	Model(t, ω, dt, N, fftp, ifftp, dispersion_term, nonlinear_function, wg.raman_model.fr, wg.γ, raman_freq_response, ω0, wg.L)
end




function _compute_error(a, b)
	sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
end

function _integrate_to_z(stepper, z, model, maxiters, reltol)
	stepper.it = 0
	while stepper.z < z

		_erk4ip_step!(stepper, model, z, reltol, maxiters)

	end
end

function simulate(u, t, model::Model; nsaves = 20, dz = 1.0, reltol = 1e-6, maxiters = 1000)

	T = t[end] - t[1]
	N = length(t)
	ν = fftshift((-N÷2:N÷2-1) ./ T)


	# initial stepper
	stepper = Stepper(model.fftp * u, model.nonlinear_function(u, model), dz, 0.0, 0.0, nothing, 0)
	zsaves = (0:nsaves) * model.L / nsaves
	dz = min(model.L / (2 * nsaves), dz)
	M = zeros(ComplexF64, (nsaves + 1, N))
	M[1, :] = stepper.U
	ϵ_hist = zeros(nsaves + 1)
	ϵ_hist[1] = 0.0


	for i ∈ 2:nsaves+1
		_integrate_to_z(stepper, zsaves[i], model, maxiters, reltol)
		ϵ_hist[i] = stepper.local_error
		M[i, :] = stepper.U
	end

	return Solution(zsaves, t, ν, ifft(M, 2), M)
end

function simulate(u, t, wg::Waveguide; args...)
	model = create_model(u, t, wg)
	simulate(u, t, model; args...)
end

function simulate(u, t, models::Vector{Model}; args...)
	sols = [simulate(u, t, models[1]; args...)]
	for model ∈ models[2:end]
		push!(sols, simulate(sols[end].At[end, :], t, model; args...))
	end
	combine(sols)
end

function simulate(u, t, wgs::Vector{Waveguide}; args...)
	models = [create_model(u, t, wg) for wg in wgs]
	simulate(u, t, models; args...)
end
