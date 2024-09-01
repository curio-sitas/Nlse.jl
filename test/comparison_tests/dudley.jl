
@testset "Supercontinuum Generation - Dudley test" begin

	n = 2^13                   # number of grid points
	twidth = 12.5e-12          # width of time window [s]
	c = 299792458              # speed of light [m/s]
	wavelength = 835e-9        # reference wavelength [m]
	w0 = (2 * pi * c) / wavelength   # reference frequency [Hz]
	dt = twidth / n
	T = (-n/2:n/2-1) .* dt # time grid

	# === input pulse
	power = 10000              # peak power of input [W]
	t0 = 28.4e-15              # duration of input [s]
	A = @. sqrt(power) * sech(T / t0) .+ 0.0im # input field [W^(1/2)]

	# === fibre parameters
	flength = 0.15             # fibre length [m]
	# betas = [beta2, beta3, ...] in units [s^2/m, s^3/m ...]
	betas = [0.0, -1.1830e-026, 8.1038e-041, -9.5205e-056, 2.0737e-070,
		-5.3943e-085, 1.3486e-099, -2.5495e-114, 3.0524e-129,
		-1.7140e-144]
	gamma = 0.11               # nonlinear coefficient [1/W/m]
	loss = 0.0                   # loss [dB/m]

	wg = Waveguide(loss, betas, gamma, wavelength, flength, raman_model = raman_linagrawaal(), self_steepening = true)
	prob = GNLSEProblem(T, wg)
	# === simulation parameters
	nsaves = 200     # number of length steps to save field at

	sol = gnlse(A, T, wg, nsaves = nsaves, dz = flength / (0.5nsaves), reltol = 1e-6)
	fn = joinpath(dirname(@__FILE__), "data\\table_dudley_test_t.csv")
	dat = CSV.read(fn, DataFrame)

	t_dudley = dat.t
	It_dudley = parse.(ComplexF64, dat.At) .|> abs2

	I = reverse(abs2.(sol.At[end, :]))

	#plot(abs.(I .- It_dudley) / maximum(I) * 100)

	err = 1 / length(I) * sum(abs.(I .- It_dudley) / maximum(I))

	@test err < 1.0 / 100

end
