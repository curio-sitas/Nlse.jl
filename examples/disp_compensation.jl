using FiberNlse
using Plots

# Simulation dimension
N = 2^13

# Fiber properties
L = 2.0e3 # Fiber length

# Signal properties
T = 100e-12 # Signal duration
λ = 1550e-9 # Wavelength
τ = 3e-12 # Pulse duration

fib1 = Waveguide(0.0, [0.0, -2.6e-26], 0.0, λ, L)
fib2 = Waveguide(0.0, [0.0, 2.6e-26], 0.0, λ, L)

t = (-N÷2:N÷2-1) * T / N



# Input construction
P₀ = 1e-3
Ψₒ = @. sqrt(P₀) / cosh(t / τ) .+ 0.0im # Soliton formula


model1 = create_model(Ψₒ, t, fib1, self_steepening = true)
model2 = create_model(Ψₒ, t, fib2)


# run the simulation
sol1 = simulate(Ψₒ, t, model1)
sol2 = simulate(sol1.At[end, :], t, model2)

sol3 = FiberNlse.combine(sol1, sol2)


heatmap(sol3.t, sol3.z, abs2.(sol3.At))
