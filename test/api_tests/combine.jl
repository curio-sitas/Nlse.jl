@testset "Combine solutions" begin

	N = 2^13

	# Fiber properties
	L = 2.0e3 # Fiber length

	# Signal properties
	T = 100e-12 # Signal duration
	λ = 1550e-9 # Wavelength
	τ = 3e-12 # Pulse duration

	fib1 = Waveguide(0.0, [0.0, -3e-26], 0.0, λ, L)
	fib2 = Waveguide(0.0, [0.0, 3e-26], 0.0, λ, L)

	t = (-N÷2:N÷2-1) * T / N



	# Input construction
	P₀ = 1e-3
	Ψₒ = @. sqrt(P₀) * sech(t / τ) .+ 0.0im # Soliton formula


	model1 = create_model(Ψₒ, t, fib1)
	model2 = create_model(Ψₒ, t, fib2)


	# run the simulation
	sol1 = simulate(Ψₒ, t, model1)
	sol2 = simulate(sol1.At[end, :], t, model2)
	sol3 = simulate(sol2.At[end, :], t, model1)
	sol4 = simulate(sol3.At[end, :], t, model2)

	sols = FiberNlse.combine([sol1, sol2, sol3, sol4])

	sola = FiberNlse.combine(sol1, sol2)
	solb = FiberNlse.combine(sol3, sol4)

	solc = FiberNlse.combine(sola, solb)


	@test solc == sols

end
