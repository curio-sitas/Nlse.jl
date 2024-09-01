

struct RamanModel
	fr::Float64
	time_response::Function
end

NoRaman = RamanModel(0.0, () -> ())

mutable struct Waveguide
	α::Union{Float64, Array{Float64}}
	βs::Array{Float64}
	γ::Union{Float64, Array{Float64}}
	λc::Float64
	L::Float64
	raman_model::RamanModel
	self_steepening::Bool
end
Waveguide(α, βs, γ, λc, L; raman_model = NoRaman, self_steepening = false) = Waveguide(α, βs, γ, λc, L, raman_model, self_steepening)


mutable struct GNLSEProblem
	ω::AbstractArray{Float64}
	dt::Float64
	N::Int
	fftp::Any
	ifftp::Any
	dispersion_term::AbstractArray{ComplexF64}
	nonlinear_function::Function
	fr::Float64
	γ::Float64
	raman_freq_response::Union{Nothing, AbstractArray{ComplexF64}}
	ω0::Float64
	L::Float64
end

mutable struct Solution
	z::AbstractArray{Float64}
	t::AbstractArray{Float64}
	f::AbstractArray{Float64}
	At::Matrix{ComplexF64}
	Af::Matrix{ComplexF64}
end

mutable struct Stepper
	U::AbstractArray{ComplexF64}
	NU::AbstractArray{ComplexF64}
	dz::Float64
	z::Float64
	local_error::Float64

	k1::Union{Nothing, AbstractArray{ComplexF64}}
	k2::Union{Nothing, AbstractArray{ComplexF64}}
	k3::Union{Nothing, AbstractArray{ComplexF64}}
	k4::Union{Nothing, AbstractArray{ComplexF64}}
	k5::Union{Nothing, AbstractArray{ComplexF64}}
	Uip::Union{Nothing, AbstractArray{ComplexF64}}
	U1::Union{Nothing, AbstractArray{ComplexF64}}
	U2::Union{Nothing, AbstractArray{ComplexF64}}
	e::Union{Nothing, AbstractArray{ComplexF64}}
	r::Union{Nothing, AbstractArray{ComplexF64}}

	it::Int
end


