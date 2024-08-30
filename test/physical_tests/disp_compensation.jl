@testset "Dispersion Compensation" begin

	# Simulation dimensionality

	N = 2^14 # Number of timesteps


	# Fiber properties

	L = 2.0e3 # Fiber length

	# Signal properties
	T = 100e-12 # Signal duration
	λ = 1550e-9 # Wavelength
	τ = 3e-12 # Pulse duration

	# Waveguides
	fib1 = Waveguide(0.0, [0.0, -2.6e-26], 0.0, λ, L)
	fib2 = Waveguide(0.0, [0.0, 2.6e-26], 0.0, λ, L)




	# Input construction
	t = (-N÷2:N÷2-1) * T / N # time cector
	P₀ = 1e-3
	Ψₒ = @. sqrt(P₀) * sech(t / τ) .+ 0.0im # Input field


	# Simulation
	sol1 = gnlse(Ψₒ, t, fib1)
	sol2 = gnlse(sol1.At[end, :], t, fib2)

	# Testing if input equals output
	@test isapprox(Ψₒ, sol2.At[end, :])

end
