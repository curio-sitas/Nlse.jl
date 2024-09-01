
"""
RamanModel 

Struct describing the raman response of the waveguide

Fields:
- `fr` is the fractional contribution of the Raman effect to the overall nonlienar effect.
- `time_response` is a user defined function returning the Raman impulse response with respect to time.

See also [raman_linagrawaal()](@ref)
"""
struct RamanModel
	fr::Float64
	time_response::Function
end


NoRaman = RamanModel(0.0, () -> ())



"""
Waveguide 

Struct describing the waveguide i.e. the propagation conditions

Fields:
- `α` is the linear optical intensity loss factor (scalar or frequency dependent)
- `βs` is the GVD orders vector (starting at 2nd order)
- `γ` is the nonlinear factor (scalar or frequency dependent)
- `λc` is the center wavelength of the problem (signal and dispersion)
- `L` is the physical length of the waveguide 
- `raman_model` is Raman model of the waveguide [RamanModel](@ref)
- `self_steepening` is a boolean to take into account the self self steepening effect

"""
mutable struct Waveguide
	α::Union{Float64, Array{Float64}}
	βs::Array{Float64}
	γ::Union{Float64, Array{Float64}}
	λc::Float64
	L::Float64
	raman_model::RamanModel
	self_steepening::Bool
end

"""
Constructs a [Waveguide](@ref) struct.

By default no Raman effect nor self steepening is activated. 
"""
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


"""
Solution

Struct containing the result of the simulation.
	
Fields: 
- `z` is the distance vector
- `t` is the duration vector 
- `f` is the frequency vector 
- `At` is the time-domain propagation field
- `Af` is the frequency-domain (spectrum) propagation field

See also [combine()](@ref)
"""
mutable struct Solution
	z::AbstractArray{Float64}
	t::AbstractArray{Float64}
	f::AbstractArray{Float64}
	At::Matrix{ComplexF64}
	Af::Matrix{ComplexF64}
end

"""
Stepper 

Struct containing the fields that are update or needed during the solver iterations.
"""
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


