

mutable struct Waveguide
	α::Any
	βs::Any
	γ::Any
	λc::Any
	L::Any
	raman_model::Any
end


Waveguide(α, βs, γ, λc, L; raman_model = :none) = Waveguide(α, βs, γ, λc, L, raman_model)



mutable struct Model
	t::Any
	ω::Any
	dt::Any
	N::Any
	fftp::Any
	ifftp::Any
	dispersion_term::Any
	nonlinear_function::Any
	fr::Any
	γ::Any
	raman_freq_response::Any
	ω0::Any
	L::Any
end

mutable struct Solution
	z::Any
	t::Any
	f::Any
	At::Any
	Af::Any
end

mutable struct Stepper
	U::Any
	NU::Any
	dz::Any
	z::Any
	local_error::Any
	k5::Any
end


struct RamanModel
	fr::Any
	time_response::Function
end
