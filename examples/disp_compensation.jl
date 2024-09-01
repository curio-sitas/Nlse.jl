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

fib1 = Waveguide(0.0, [0.0, -3e-26], 0.0, λ, L)
fib2 = Waveguide(0.0, [0.0, 3e-26], 0.0, λ, L)

t = (-N÷2:N÷2-1) * T / N



# Input construction
P₀ = 1e-3
Ψₒ = @. sqrt(P₀) * sech(t / τ) .+ 0.0im # Soliton formula


model1 = GNLSEProblem(t, fib1)
model2 = GNLSEProblem(t, fib2)


# run the simulation
sol1 = gnlse(Ψₒ, t, fib1)
sol2 = gnlse(sol1.At[end, :], t, fib1)
sol3 = gnlse(sol2.At[end, :], t, model2)
sol4 = gnlse(sol3.At[end, :], t, model2)

sols = FiberNlse.combine([sol1, sol2, sol3, sol4])

sols2 = gnlse(Ψₒ, t, [fib1, fib2, fib1, fib2])

heatmap(sols2.t, sols2.z, (abs2.(sols2.At)))
