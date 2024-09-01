@testset "NLSE Reversibility" begin
	# Simulation dimension
	N = 2^13

	# Fiber properties
	L = 2.0e3 # Fiber length

	# Signal properties
	T = 100e-12 # Signal duration
	λ = 1550e-9 # Wavelength
	τ = 15e-12 # Pulse duration

	fib1 = Waveguide(1e-5, [0.0, -2.6e-26], 1.1e-3, λ, L)
	fib2 = Waveguide(-1e-5, [0.0, 2.6e-26], -1.1e-3, λ, L)

	t = (-N÷2:N÷2-1) * T / N



	# Input construction
	P₀ = 0.5
	Ψₒ = @. sqrt(P₀) * sech(t / τ) .+ 0.0im # Soliton formula


	# run the simulation
	sol1 = gnlse(Ψₒ, t, fib1)
	sol2 = gnlse(sol1.At[end, :], t, fib2)

	# Testing signal propagation (including losses)
	@test FiberNlse._compute_error(sol1.At[1, :] .|> abs2, sol2.At[end, :] .|> abs2) < 1 / 100

end
