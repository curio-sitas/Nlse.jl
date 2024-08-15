using FFTW
using StatsBase: mean
using LinearAlgebra
c = 3e8
mutable struct Material
	L::Real
	α::Real
	β::Vector{Real}
	γ::Real
	λ::Real
	fr::Real
	τl::Real
	τvib::Real
end

smf28(L::Real, λ::Real) = Material(L, 0.046e-3, [0.0, -2.1682619391414893e-26], 1.1e-3, λ, 0.18, 32e-15, 12.2e-15)

struct NLSEProblem
	t::Vector{Real}
	u0::Vector{ComplexF64}
	material::Material
	##
end

function hr(material, t)
	@. 0.0 + (t > 0) * (1.0 / material.τvib^2 + 1.0 / material.τl^2) * material.τvib * exp(-t / material.τl) * sin(t / material.τvib)
end

function solve(prob::NLSEProblem, dz::Real; n_saves = nothing, tol::Real = 1e-3)

	N = length(prob.t)
	T = prob.t[end] - prob.t[1]

	ν = fftshift(collect(-N/2:N/2-1)) / T
	global ν0 = c / prob.material.λ

	L = prob.material.L
	α = prob.material.α
	γ = prob.material.γ

	function f_NLSE(u, U)
		IT = u .* abs2.(u)
		1im * γ * fft(IT .+ 1im * ifft(1im * 2pi * ν .* fft(IT)) ./ (2pi * ν0))
	end

	d = -0.5α * ones(N)

	for i in eachindex(prob.material.β)
		d = d .+ 1im * prob.material.β[i] * (2π * ν) .^ (i) / factorial(i)
	end

	u = ComplexF64.(prob.u0)
	U = fft(u)
	NU = f_NLSE(u, U)
	zk = 0.0


	if !isnothing(n_saves)
		z_saves = LinRange(0, L, n_saves)
		us = zeros(ComplexF64, n_saves, N)
	else
		z_saves = [L]
		us = zeros(ComplexF64, 1, N)
	end



	for (i, zi) in enumerate(z_saves)
		while zk < zi

			e = exp.(0.5 * dz * d)
			Uip = e .* U

			k1 = e .* NU
			U2 = Uip .+ 0.5 * dz * k1
			u2 = ifft(U2)
			k2 = f_NLSE(u2, U2)

			U3 = Uip .+ 0.5 * dz * k2
			u3 = ifft(U3)
			k3 = f_NLSE(u3, U3)

			U4 = e .* (Uip .+ dz * k3)
			u4 = ifft(U4)
			k4 = f_NLSE(u4, U4)

			r = e .* (Uip .+ dz * (k1 / 6.0 .+ k2 / 3.0 .+ k3 / 3.0))
			U1 = r .+ dz * k4 / 6.0
			u1 = ifft(U1)

			k5 = f_NLSE(u1, U1)

			U2 = r .+ dz * (k4 / 15.0 .+ k5 / 10.0)

			err = _compute_error(U1, U2)

			dzopt = max(0.5, min(2.0, 0.9 * sqrt(sqrt(tol / err)))) * dz

			if err <= tol

				dz = min(dzopt, abs(zk - zi))
				zk = zk + dz

				U = U1
				u = u1
				NU = k5

			else
				dz = dzopt
			end
		end


		us[i, :] = u

	end


	return z_saves, us

end

_compute_error(a, b) = norm(a - b) / norm(a)

function dBm2W(PdBm)
	10^(0.1 * PdBm - 3)
end

function W2dBm(Pw)
	10 * log10(Pw) + 30
end

output(sol::Matrix{ComplexF64}) = sol[:, end]

gaussian(t, Δt, C) = @. exp(-2 * log(2) * (1.0 + 1im * C) * (t / Δt)^2)


function fwhm(t, u)

	idx = findall(L1normalize(u) .>= 0.5)
	return abs(t[idx[end]] - t[idx[1]])


end

L1normalize(x) = x ./ maximum(x)
L2normalize(x) = x ./ sqrt(sum(abs2.(x)))
RMSnormalize(x) = x ./ sqrt(mean(abs2.(x)))

