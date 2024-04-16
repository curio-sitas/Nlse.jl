using FiberNlse

T = 100e-12
N = 2^13
L = 1000.0
λ = 1550e-9
P0 = 1.0
τ = 5.0e-12

t = T * (-N÷2:N÷2-1) / N

u0 = @. sqrt(P0) * sech(t / τ) + 0im

SMFmat = Material([0.0, 0.0], 0.0, 0.0)
smf = Fiber(SMFmat, L, λ)

sol = propagate(u0, T, smf)

