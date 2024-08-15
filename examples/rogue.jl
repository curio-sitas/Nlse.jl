using Revise


using FiberNlse, Plots, StatsBase, DSP

FiberNlse.setPhaseConvention(:positive)

D = 18e-6
α = 0*0.026e-3
γ = 1.1e-3
L = 50000

# Signal properties
λ = 1550e-9 # Wavelength

Nt, Nz = (3000,100)
f0 = 10e9
T = 3/f0
t0 = 500e-12
t = LinRange(-1,1,Nt)*T/2
C = -0.8
P0 = 0.3
mo=0.01
σ=0*1e-4
n= σ*randn(Nt)
S = sqrt(P0).*(1.0.+mo*cos.(2pi*f0*t) .+ n)


fib = Fiber(L, dispersion(D, λ), γ, α, λ)
field = propagate(S, fib, T, Nz)

begin
    plot(t/ps,field.ψ[end, :] .|> abs2, label="end")
    plot!(t/ps,S .|> abs2, label="begin")

end

heatmap(abs2.(field.ψ))

@gif for i in 1:Nz
f,sp = FiberNlse.spectrum(field.ψ[i,:],T/Nt)
plot(f/1e9,10log10.(sp .|> abs2),ylims=(-200,0),label="$(i)/$(Nz)")
end