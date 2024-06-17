using Revise
using Plots
using DSP: unwrap

using FFTW
using ProgressBars
using StatsBase: mean


mutable struct Material
	L::Real
	α::Real
	β::Vector{Real}
	γ::Real
	λ::Real
end

mutable struct Model
	D::Function
	N::Function
	FFT::Any
	IFFT::Any
end

function Model()
	Model((x) -> nothing, (x) -> nothing, nothing, nothing)
end

function construct_model(mat::Material, model::Symbol = :NLSE)
	@info "Not implemented"
end

smf28(L::Real, λ::Real) = Material(L, 0.046e-3, [0.0, -2.1682619391414893e-26], 1.1e-3, λ)

struct NLSEProblem
	t::Vector{Real}
	u0::Vector{ComplexF64}
	material::Material
	##
end

function _compute_error(a, b)
	sqrt(sum(abs2.(a .- b)) ./ sum(abs2.(a)))
end

function _erk4ip_step!(û, k5, z, dz, model::Model)

	e = exp.(0.5 * dz * model.D(z))
	Uip = e .* û

	if isnothing(k5)
		k1 = e .* model.N(model.IFFT * û, model.FFT)
	else
		k1 = model.D(z) .* k5
	end

	k2 = model.N(model.IFFT * (Uip .+ 0.5 * dz * k1), model.FFT)

	k3 = model.N(model.IFFT * (Uip .+ 0.5 * dz * k2), model.FFT)

	k4 = model.N(model.IFFT * (e .* (Uip .+ dz * k3)), model.FFT)

	r = e .* (Uip .+ dz * (k1 / 6.0 .+ k2 / 3.0 .+ k3 / 3.0))

	U1 = r .+ dz * k4 / 6.0

	k5 = model.N(model.IFFT * U1, model.FFT)

	U2 = r .+ dz * (k4 / 15.0 .+ k5 / 10.0)

	return U1, _compute_error(U1, U2) #MSE

end

function solve(prob::NLSEProblem, dz::Real; tol::Real = 1e-6, fullsol::Bool = true, showprogress::Bool = false, Ladim::Real = 0.0)

	# Define constants and vectors
	N = length(prob.t)
	T = prob.t[end] - prob.t[1]
	ν = fftshift(collect(-N/2:N/2-1))
	dz = dz / prob.material.L
	Pmax = maximum(abs2.(prob.u0))
	α = prob.material.α * prob.material.L
	γ = prob.material.γ * prob.material.L * Pmax
	u0 = prob.u0 ./ sqrt(Pmax)

	D = -0.5α * ones(N)
	for i in eachindex(prob.material.β)
		D = D .+ 1im * prob.material.β[i] * (2π * ν) .^ (i) / factorial(i) / T^i * prob.material.L
	end

	#!TODO fft flags
	#!TODO construct_model for NLSE and GNLSE



	NL(u, FFT) = 1im * γ * (FFT * (u .* abs2.(u)))
	model = Model((x) -> D, (u, v) -> NL(u, v), FFTW.plan_fft(prob.u0), FFTW.plan_ifft(prob.u0))


	z = [0.0]


	U = fft(u0)
	M = U
	zk = 0
	k5 = nothing
	ϵ_hist = []

	stepper = Stepper(fft(u0), nothing, zk, dz, 0.0)

	while zk < 1.0

		_erk4ip_step!(U, k5, zk, dz, model)

		dzopt = max(0.5, min(2.0, 0.9 * sqrt(sqrt(tol / err)))) * dz

		if err <= tol

			zk = zk + dz
			stepper.z = zk
			stepper.dz = min(dzopt, L - zk)

			U = U1
			append!(z, zk)
			append!(M, U)
			append!(ϵ_hist, stepper.ϵ)


		else
			dz = dzopt
		end
	end

	return errs, z * prob.material.L, ifft(reshape(M, N, length(M) ÷ N), 1) .* sqrt(Pmax)


end


output(sol::Matrix{ComplexF64}) = sol[:, end]

gaussianPulse(t, Δt) = exp.(-4 * log(2) * (t / Δt) .^ 2)

begin

	ps = 1e-12
	nm = 1e-9
	km = 1e3


	N = 2^15
	T = 100ps
	t = T * (-N÷2:N÷2-1) / N

	λ = 1550nm
	L = 5e3


	fiber = smf28(L, λ)
	fiber.α = 0 * 0.000045

	Δt = 19 * ps
	P0 = abs(fiber.β[2]) / (fiber.γ * Δt^2)
	u0 = sqrt(P0) * sech.(t / Δt)

	prob = NLSEProblem(t, u0, fiber)
	err, z, sol = solve(prob, 1.0, tol = 1e-6)

	heatmap(t / ps, z / km, abs2.(sol'), levels = 5)
end

plot(t / ps, abs2.(u0))
plot!(t / ps, abs2.(output(sol)))

plot(cumsum(err))
hline!([1e-6])




mutable struct Stepper
	u::Vector{ComplexF64}
	k5::Vector{ComplexF64}
	z::Float64
	dz::Float64
	ϵ::Float64
end

