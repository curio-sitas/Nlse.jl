using FiberNlse
using Plots


ps = 1e-12
km = 1e3

β2 = -20 * ps^2 / km
α = 0.046 / km
γ = 1.1 / km

λ = 1.55e-6
L = 10e3

fiber = Waveguide(α, [0.0, β2], γ, λ, L, raman_linagrawaal())


# Signal 
N = 2^14
T = 100ps
t = T * (-N÷2:N÷2-1) / N

# Params

P0 = 0.5
Vpirf = 4.2
Vpip = 4.0
Vpidc = 6.6
Δθ = pi
Ω0 = 2pi * 10e9
Vb = 3.3
V0 = 2.1
q = 1.5


U = @. sqrt(P0) * cos(0.5pi * V0 / Vpirf * cos(Ω0 * t) + 0.5pi * Vb / Vpidc) * cis(pi * q * V0 / Vpip * cos(Ω0 * t + Δθ))
model = create_model(U, t, fiber; raman = true, self_steepening = true)

sol = simulate(U, t, model, nsaves = 200)
using FFTW
heatmap(abs2.(fftshift(sol.At, 2)))

zindex = findall(sol.z .% 1000 .== 0)

plot(camera = (50, 20))
for zi ∈ zindex
	plot!(t * 1e12, 1e-3 * sol.z[zi] * ones(N), fftshift(sol.At[zi, :]) .|> abs2, color = :black, lw = 2, label = false)
end
zvals = sol.z[zindex] / 1000
plot!(yticks = zvals, ylims = (0, zvals[end]))
plot!(xticks = (-50:19:50), xlims = (-50, 50))
display(current())

surface(abs2.(fftshift(sol.At, 2)))
