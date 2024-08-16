@testset "Self-Phase Modulation" begin

	# Simulation dimension
	N = 2^14

	# Fiber properties
	L = 5.0e3 # Fiber length

	# Signal properties
	τ = 10e-12 # Pulse duration
	T = 100e-12 # Signal duration
	λ = 1550e-9 # Wavelength
	n = 1 # Soliton number

	fib = Waveguide(0.0, [0.0, 0.0], 1.1e-3, λ, L)

	t = (-N÷2:N÷2-1) * T / N


	# Input construction
	P₀ = 1.0 # Soliton power
	Ψₒ = sqrt(P₀) * cos.(t ./ τ) .+ 0.0im# Soliton formula

	model = create_model(Ψₒ, t, fib)

	sol = simulate(Ψₒ, t, model, nsaves = 2^5, reltol = 1e-9, dz = 1e-6 * L)

	origin(x) = x .- x[1]

	φ = unwrap(angle.(sol.At[end, :]), range = pi) |> origin
	φₜₕ = (fib.γ * abs2.(Ψₒ) * L) |> origin


	@test FiberNlse._compute_error(φₜₕ, φ) < 1 / 100

end
