

struct RamanModel
	fr::Float64
	time_response::Function
end

NoRaman = RamanModel(0.0, () -> ())

mutable struct Waveguide{T <: Real}
	α::Union{T, Array{T}}
	βs::Array{T}
	γ::Union{T, Array{T}}
	λc::T
	L::T
	raman_model::RamanModel
end
Waveguide(α, βs, γ, λc, L; raman_model = NoRaman) = Waveguide(α, βs, γ, λc, L, raman_model)


mutable struct Model{T <: Real}
	t::AbstractArray{T}
	ω::AbstractArray{T}
	dt::T
	N::Int
	fftp::Any
	ifftp::Any
	dispersion_term::AbstractArray{ComplexF64}
	nonlinear_function::Function
	fr::T
	γ::T
	raman_freq_response::Union{Nothing, AbstractArray{ComplexF64}}
	ω0::T
	L::T
end

mutable struct Solution{T <: Real}
	z::AbstractArray{T}
	t::AbstractArray{T}
	f::AbstractArray{T}
	At::Matrix{ComplexF64}
	Af::Matrix{ComplexF64}
end

mutable struct Stepper{T <: Real}
	U::AbstractArray{ComplexF64}
	NU::AbstractArray{ComplexF64}
	dz::T
	z::T
	local_error::T
	k5::Union{Nothing, AbstractArray{ComplexF64}}
end


