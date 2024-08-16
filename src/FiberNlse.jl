module FiberNlse

using FFTW

export simulate, Waveguide, raman_linagrawaal, create_model, Solution

export combine

c = 299792458

@info "The integrated signal is supposed to follow the negative phase convention."

include("datatypes.jl")
include("api.jl")
include("algorithm.jl")
include("utils.jl")

end

