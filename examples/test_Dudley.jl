using FiberNlse
using Plots
using FFTW
n = 2^13;                   # number of grid points
twidth = 12.5e-12;          # width of time window [s]
c = 299792458;              # speed of light [m/s]
wavelength = 835e-9;        # reference wavelength [m]
w0 = (2 * pi * c) / wavelength;   # reference frequency [Hz]
dt = twidth / n;
T = (-n/2:n/2-1) .* dt; # time grid

# === input pulse
power = 10000;              # peak power of input [W]
t0 = 28.4e-15;              # duration of input [s]
A = @. sqrt(power) * sech(T / t0) .+ 0.0im; # input field [W^(1/2)]

# === fibre parameters
flength = 0.15;             # fibre length [m]
# betas = [beta2, beta3, ...] in units [s^2/m, s^3/m ...]
betas = [0.0, -1.1830e-026, 8.1038e-041, -9.5205e-056, 2.0737e-070,
	-5.3943e-085, 1.3486e-099, -2.5495e-114, 3.0524e-129,
	-1.7140e-144];
gamma = 0.11;               # nonlinear coefficient [1/W/m]
loss = 0.0;                   # loss [dB/m]

wg = Waveguide(loss, betas, gamma, wavelength, flength, raman_model = raman_linagrawaal(), self_steepening = true)
prob = GNLSEProblem(T, wg);
# === simulation parameters
nsaves = 200;     # number of length steps to save field at

@time sol = gnlse(A, T, wg, nsaves = 200, dz = flength / (2 * nsaves), reltol = 1e-5);

# Spectrum

WL = c ./ (reverse(fftshift(sol.f)) .+ c ./ wavelength);
iis = findall(450e-9 .< WL .< 1350e-9); # wavelength grid
lIW = 10 * log10.(abs.(reverse(ifftshift(sol.Af, 2), dims = 2) .^ 2)); # log scale spectral intensity
mlIW = lIW .- maximum(lIW);       # max value, for scaling plot

heatmap(WL[iis] * 1e9, sol.z, mlIW[:, iis], clims = (-40, 0), xlims = (450, 1350))
xlabel!("Wavelength / nm")
ylabel!("Distance / m")

# Timeseries

lIT = 10 * log10.(abs.(reverse(sol.At, dims = 2)) .^ 2); # log scale temporal intensity
mlIT = lIT .- maximum(lIT);       # max value, for scaling plot
heatmap(T .* 1e12, sol.z, mlIT, clims = (-40, 0), xlims = (-0.5, 5))   # plot as pseudocolor map

xlabel!("Delay / ps");
ylabel!("Distance / m");
