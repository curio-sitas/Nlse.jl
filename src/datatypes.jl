

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
end
Waveguide(α, βs, γ, λc, L; raman_model = NoRaman) = Waveguide(α, βs, γ, λc, L, raman_model)


mutable struct Model
	t::AbstractArray{Float64}
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
	k5::Union{Nothing, AbstractArray{ComplexF64}}
end


