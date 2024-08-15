@testset "Self-Phase Modulation" begin

	#TODO : rewrite clean test

	# Simulation dimension
	Nₜ = 2^13

	# Fiber properties
	L = 5.0e3 # Fiber length

	# Signal properties
	τ = 40e-12 # Pulse duration
	T = 4000e-12 # Signal duration
	λ = 1550e-9 # Wavelength
	N = 1 # Soliton number

	α = 0.0

	fib = Waveguide(α, [0.0, -2.6e-26], 1.1e-3, λ, L)

	t = (-Nₜ÷2:Nₜ÷2-1) * T / Nₜ


	# Input construction
	P₀ = abs((fib.βs[2] / fib.γ / τ^2) * N^2) # Soliton power
	Ψₒ = sqrt(P₀) * sech.(t ./ τ) .+ 0.0im# Soliton formula

	model = create_model(Ψₒ, t, fib)

	sol = simulate(Ψₒ, t, model, nsaves = 200)

	φ = unwrap(angle.(sol.At[end, :]), range = pi)
	φₜₕ = -fib.γ * abs2.(Ψₒ) * L

	@test (sum(abs.(φₜₕ .- φ) ./ φₜₕ) / length(φₜₕ)) < 0.001
end
