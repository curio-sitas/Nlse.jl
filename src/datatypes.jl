struct Material

	β::Vector{Real}
	α::Real
	γ::Real

end

struct Fiber

	mat::Material
	L::Real
	λ::Real

end

struct PropagField

	u::AbstractVecOrMat{ComplexF64}
	û::AbstractVecOrMat{ComplexF64}
	z::Vector{Real}
	t::Vector{Real}

end
