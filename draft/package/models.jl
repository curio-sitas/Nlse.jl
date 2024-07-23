


mutable struct Model
    FFT::Any
    IFFT::Any
    dt::Any
    fr::Any
    γ::Any
    RW::Any
    ω::Any
    ω0::Any
    D::Any
    NL::Any
end

mutable struct Stepper
    U::Any
    NU::Any
    dz::Any
    z::Any
    local_error::Any
    k5::Any
end
