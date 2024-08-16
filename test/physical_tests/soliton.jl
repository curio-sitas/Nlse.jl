@testset "Soliton Invariance" begin

	# Simulation dimension
	N = 2^13

	# Fiber properties
	L = 5.0e3 # Fiber length

	# Signal properties
	τ = 40e-12 # Pulse duration
	T = 4000e-12 # Signal duration
	λ = 1550e-9 # Wavelength
	n = 1 # Soliton number

	α = 0.0

	fib = Waveguide(α, [0.0, -2.6e-26], 1.1e-3, λ, L)

	t = (-N÷2:N÷2-1) * T / N


	# Input construction
	P₀ = abs((fib.βs[2] / fib.γ / τ^2) * n^2) # Soliton power
	Ψₒ = sqrt(P₀) * sech.(t ./ τ) .+ 0.0im # Soliton formula

	model = create_model(Ψₒ, t, fib)

	sol = simulate(Ψₒ, t, model, nsaves = 200)

	# Testing soliton propagation (including losses)
	@test isapprox(abs2.(Ψₒ .* exp(-0.5 * fib.α * L)), abs2.(sol.At[end, :]), atol = 1e-4)
end
