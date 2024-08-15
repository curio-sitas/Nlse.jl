# IMPORTS 
using Revise
using FFTW, Plots
using ProfileView

# MAIN
begin
	# CONSTANTS
	c = 299792458

	# WAVEGUIDE
	fr = 0.18
	τl = 32e-15
	τvib = 12.2e-15
	L = 0.5
	α = 0.0
	γ = 0.11
	λc = 835e-9
	ω0 = 2π * c / λc
	βs = [0.0, -1.1830e-026]

	# SIGNAL 
	T = 12.5e-12 # Duration
	N = 2^14 # Time samples
	τ = 0.05e-12 / (2log(1 + sqrt(2)))
	P = abs(βs[2] / (γ * τ^2))
	n = 3
	# VECTORS 
	dt = T / N
	t = T * (-N÷2:N÷2-1) / N
	u0 = n * sqrt(P) * sech.(t ./ τ)
end









wg = Waveguide(α, βs, γ, λc, L, raman_linagrawaal())
model = create_model(u0, t, wg; self_steepening = true, raman = true, normalize = false);


SOL = simulate(u0, t, model; dz = 0.01, nsaves = 256, reltol = 1e-6, maxiters = 100)


TT = -1e-12 .<= t .<= 2e-12
SOLr = reverse(SOL.At, dims = 2)[:, TT] .|> abs2
SOLr = SOLr ./ maximum(SOLr)
MIT = 10log10.(SOLr)

heatmap(t[TT] * 1e12, SOL.z, MIT, clims = (-40, 0))
xlabel!("Time [ps]")
ylabel!("Distance [m]")
