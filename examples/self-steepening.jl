using Revise
using Plots
using  FiberNlse
# Simulation dimension
Nₜ, Nₗ = (1000, 1000);

# Fiber properties
L = 5e3; # Fiber length

# Signal properties
λ = 1550e-9; # Wavelength
τ = 20e-12; # Pulse duration
T = 20*τ; # Signal duration
N = 2# Soliton number

fib = FiberNlse.Fiber(
    L,
    FiberNlse.Dispersion([-2.6e-26]),
    1e-3,
    1e-5,
    λ,
);
fib.D
t = (0:(Nₜ-1)) * T / Nₜ .- 0.5T

# Input construction
P₀ =1#abs((fib.D.β[1]/fib.γ/τ^2)*N^2) # Soliton power
Ψₒ = @. sqrt(P₀) * exp(-0.5(t / τ)^2) # Soliton formula

field = FiberNlse.propagate(Ψₒ, [fib], T, Nₗ; progress = true); # run the simulation

begin
    plot(t*1e12,abs2.(field.ψ[end, :]))
    plot!(t*1e12,abs2.(Ψₒ))
    xlims!(10τ*1e12.*(-1,1))
end
sum(abs2.(Ψₒ))
sum(abs2.(field.ψ[end, :]))

# Testing soliton propagation (including losses)
@test isapprox(abs2.(Ψₒ .* exp(-fib.α * L)), abs2.(field.ψ[end, :]), atol = 1e-5)

heatmap(abs2.(field.ψ))