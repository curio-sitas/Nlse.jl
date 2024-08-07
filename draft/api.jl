# IMPORTS 
using Revise
using FFTW, Plots


# MAIN
begin
	# CONSTANTS
	c = 299792458

	# WAVEGUIDE
	fr = 0.18
	τl = 32e-15
	τvib = 12.2e-15
	L = 0.5
	α = 0.0
	γ = 0.11
	λc = 835e-9
	ω0 = 2π * c / λc
	βs = [0.0, -1.1830e-026]

	# SIGNAL 
	T = 12.5e-12 # Duration
	N = 2^14 # Time samples
	τ = 0.05e-12 / (2log(1 + sqrt(2)))
	P = abs(βs[2] / (γ * τ^2))
	n = 3
	# VECTORS 
	dt = T / N
	t = T * (-N÷2:N÷2-1) / N
	u0 = n * sqrt(P) * sech.(t ./ τ)
end


ν = fftshift((-N÷2:N÷2-1) ./ T);
ω = 2π .* ν;


function NL_spm(u, model)
	return 1.0im .* model.γ * model.fftp * (u .* abs2.(u))
end

function NL_spm_self_steepening(u, model)
	η = (1 .+ model.ω / model.ω0)
	return 1.0im .* model.γ * η .* (model.fftp * (u .* abs2.(u)))
end

function NL_spm_self_steepening_raman(u, model)
	η = (1 .+ model.ω / model.ω0)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	M = model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS))
	return 1.0im .* model.γ * η .* M
end

function NL_spm_raman(u, model)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	M = model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS))
	return 1.0im .* model.γ .* M
end

function setup_nonlinearity_term(self_steepening::Bool = true, raman::Bool = true)

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


struct RamanModel
	fr::Any
	time_response::Function
end


function raman_linagrawaal(fr = 0.18, τl = 32e-15,
	τvib = 12.2e-15)

	function time_response(t)
		hr = 0.0 .+ (t .>= 0) .* ((τvib^2 + τl^2) / τvib / τl^2 * exp.(-t / τl) .* sin.(t / τvib))
		return hr
	end

	return RamanModel(fr, time_response)

end

mutable struct Waveguide
	α::Any
	βs::Any
	γ::Any
	λc::Any
	L::Any
	raman_model::Any
end
Waveguide(α, βs, γ, λc, L; raman_model = :none) = Waveguide(α, βs, γ, λc, L, raman_model)


function create_model(u, t, wg::Waveguide; self_steepening::Bool = false, raman::Bool = false, normalize = true)

	fftp = plan_fft(ComplexF64.(u))
	ifftp = plan_ifft(ComplexF64.(u))

	T = t[end] - t[1]
	N = length(u)
	dt = T / N

	ν = fftshift((-N÷2:N÷2-1) ./ T)
	ω = 2π .* ν

	ω0 = 2pi * c / wg.λc


	raman_freq_response = nothing
	if raman && wg.raman_model != :none
		raman_freq_response = conj((fftp * ifftshift(wg.raman_model.time_response(t))))
	end

	dispersion_term = -0.5α .+ 1im * sum([(wg.βs[i] / factorial(i)) .* ω .^ i for i in eachindex(wg.βs)])
	nonlinear_function = setup_nonlinearity_term()
	model = Model(t, ω, dt, N, fftp, ifftp, dispersion_term, nonlinear_function, wg.raman_model.fr, wg.γ, raman_freq_response, ω0, wg.L)
end

mutable struct Model
	t::Any
	ω::Any
	dt::Any
	N::Any
	fftp::Any
	ifftp::Any
	dispersion_term::Any
	nonlinear_function::Any
	fr::Any
	γ::Any
	raman_freq_response::Any
	ω0::Any
	L::Any
end

function _compute_error(a, b)
	sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
end

function _integrate_to_z(stepper, z, model, maxiters, reltol)
	it = 0
	while stepper.z < z

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
			it = it + 1
			if (it >= maxiters)
				throw(ErrorException("Max number of iteration exceeded!"))
			end
		end
	end
end

mutable struct Solution
	z::Any
	t::Any
	f::Any
	At::Any
	Af::Any
end

mutable struct Stepper
	U::Any
	NU::Any
	dz::Any
	z::Any
	local_error::Any
	k5::Any
end


function simulate(u, t, model; nsaves = 20, dz = 1.0, reltol = 1e-6, maxiters = 1000)

	T = t[end] - t[1]
	ν = fftshift((-N÷2:N÷2-1) ./ T)


	# initial stepper
	stepper = Stepper(model.fftp * u, model.nonlinear_function(u, model), dz, 0.0, 0.0, 0.0)
	zsaves = (0:nsaves) * L / nsaves
	dz = min(L / (2 * nsaves), dz)
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



wg = Waveguide(α, βs, γ, λc, L, raman_linagrawaal())
model = create_model(u0, t, wg; self_steepening = true, raman = true, normalize = false)


SOL = simulate(u0, t, model; dz = L / 100, nsaves = 256, reltol = 1e-6, maxiters = 100)


TT = -1e-12 .<= t .<= 2e-12
SOLr = reverse(SOL.At, dims = 2)[:, TT] .|> abs2
SOLr = SOLr ./ maximum(SOLr)
MIT = 10log10.(SOLr)

heatmap(t[TT] * 1e12, SOL.z, MIT, clims = (-40, 0))
xlabel!("Time [ps]")
ylabel!("Distance [m]")
