using Revise


using FiberNlse, Plots, StatsBase, DSP

FiberNlse.setPhaseConvention(:positive)

D = 18e-6
α = 0.026e-3
γ = 1.1e-3
L = 1000.0

# Signal properties
λ = 1550e-9 # Wavelength

Nt, Nz = (5000,200)
T = 300e-12
t0 = 20e-12
t = LinRange(-1,1,Nt)*T/2
C = -0.8
P0 = 0.2
y = exp.(-((t./t0).^2) .*(1-im*C))
S = sqrt(P0)*y./sqrt(mean(abs2.(y)))

fib = Fiber(L, dispersion(D, λ), γ, α, λ)
field = propagate(S, fib, T, Nt)

φ = instFreq(field.ψ[end, :],t)/1e9
φₜₕ = instFreq(S, t)/1e9

begin
    plot(t/ps,φ, label="end")
    plot!(t/ps,φₜₕ, label="begin")
    xlims!(-50,50)
    ylims!(-20,20)
end

begin
    plot(t/ps,field.ψ[end, :] .|> abs2, label="end")
    plot!(t/ps,S .|> abs2, label="begin")
    xlims!(-50,50)
end


S = sqrt(P0)*(noise+cos(2pi*f0*t))