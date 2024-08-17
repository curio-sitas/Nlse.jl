module FiberNlse

using FFTW

export Waveguide, Solution, Model, RamanModel

export simulate, raman_linagrawaal, create_model

export combine

# Light celerity
c = 299792458

@info "The integrated signal is supposed to follow the negative phase convention."

include("datatypes.jl")
include("api.jl")
include("raman.jl")
include("nonlinearity.jl")
include("algorithm.jl")
include("utils.jl")

end

