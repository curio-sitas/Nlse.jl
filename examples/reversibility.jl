using FiberNlse
using Plots
# Simulation dimension
N = 2^15

# Fiber properties
L = 6.0e3 # Fiber length

# Signal properties
T = 500e-12 # Signal duration
λ = 1550e-9 # Wavelength
τ = 10e-12 # Pulse duration

fib1 = Waveguide(0.0, [0.0, -2.6e-26], 1.1e-3, λ, L)
fib2 = Waveguide(0.0, [0.0, 2.6e-26], -1.1e-3, λ, L)

t = (-N÷2:N÷2-1) * T / N



# Input construction
P₀ = 0.8
Ψₒ = @. sqrt(P₀) * sech(t / τ) .+ 0.0im # Soliton formula

model1 = GNLSEProblem(t, fib1)
model2 = GNLSEProblem(t, fib2)


# run the simulation
sol1 = gnlse(Ψₒ, t, model1)
sol2 = gnlse((sol1.At[end, :]), t, model2)

sol3 = combine(sol1, sol2)

heatmap(abs2.(sol3.At))

plot(sol1.At[1, :] .|> abs2)
plot!(sol2.At[end, :] .|> abs2)
