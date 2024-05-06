using FFTW
using ProgressBars
using StatsBase: mean
PB_STEPS = 1000

mutable struct Material
    L::Real
    α::Real
    β::Vector{Real}
    γ::Real
    λ::Real
end

smf28(L::Real, λ::Real) = Material(L, 0.046e-3, [0.0, -2.1682619391414893e-26], 1.1e-3, λ)

struct NLSEProblem
    t::Vector{Real}
    u0::Vector{Complex}
    material::Material
    ##
end

function solve(prob::NLSEProblem, dz::Real; tol::Real=1e-6, fullsol::Bool=true, showprogress::Bool=false, Ladim::Real=0.0)

    N = length(prob.t)
    T = prob.t[end] - prob.t[1]

    ν = fftshift(collect(-N/2:N/2-1))
    t = prob.t / T

    if Ladim == 0.0
        Ladim = prob.material.L
    end

    L = prob.material.L / Ladim
    Pmax = maximum(abs2.(prob.u0))
    dz = dz / Ladim
    α = prob.material.α * Ladim
    γ = prob.material.γ * Ladim * Pmax

    function f_NLSE(u, U)
        1im * γ * fft(u .* abs2.(u))
    end

    z = [0.0]

    u0 = prob.u0 ./ sqrt(Pmax)
    d = -0.5α * ones(N)

    if showprogress
        progress = ProgressBar(total=PB_STEPS, printing_delay=0.01)
    end

    for i in eachindex(prob.material.β)
        d = d .+ 1im * prob.material.β[i] * (2π * ν) .^ (i) / factorial(i) / T^i * Ladim
    end

    u = u0
    U = fft(u)
    NU = f_NLSE(u, U)
    M = u0
    zk = 0
    c = 0
    while zk < L

        e = exp.(0.5 * dz * d)
        Uip = e .* U

        k1 = e .* NU
        U2 = Uip .+ 0.5 * dz * k1
        u2 = ifft(U2)
        k2 = f_NLSE(u2, U2)

        U3 = Uip .+ 0.5 * dz * k2
        u3 = ifft(U3)
        k3 = f_NLSE(u3, U3)

        U4 = e .* (Uip .+ dz * k3)
        u4 = ifft(U4)
        k4 = f_NLSE(u4, U4)

        r = e .* (Uip .+ dz * (k1 / 6.0 .+ k2 / 3.0 .+ k3 / 3.0))
        U1 = r .+ dz * k4 / 6.0
        u1 = ifft(U1)

        k5 = f_NLSE(u1, U1)

        U2 = r .+ dz * (k4 / 15.0 .+ k5 / 10.0)

        err = sqrt(sum(abs2.(U1 .- U2)) ./ sum(abs2.(U1))) #MSE

        dzopt = max(0.5, min(2.0, 0.9 * sqrt(sqrt(tol / err)))) * dz

        if err <= tol

            if showprogress
                update(progress, round(Int, PB_STEPS * dz))
            end

            zk = zk + dz
            z = vcat(z, zk)
            dz = min(dzopt, L - zk)

            U = U1
            u = u1
            NU = k5

            if fullsol
                append!(M, u)
            end



        else
            dz = dzopt
        end
    end



    if fullsol == true
        return z * Ladim, reshape(M, N, length(M) ÷ N) .* sqrt(Pmax)
    else
        return z * Ladim, u .* sqrt(Pmax)
    end

end


function dBm2W(PdBm)
    10^(0.1 * PdBm - 3)
end

function W2dBm(Pw)
    10 * log10(Pw) + 30
end

output(sol::Matrix{ComplexF64}) = sol[:, end]

gaussian(t, Δt, C) = @. exp(-2 * log(2) * (1.0 + 1im * C) * (t / Δt)^2)


function fwhm(t, u)

    idx = findall(L1normalize(u) .>= 0.5)
    return abs(t[idx[end]] - t[idx[1]])


end

L1normalize(x) = x ./ maximum(x)
L2normalize(x) = x ./ sqrt(sum(abs2.(x)))
RMSnormalize(x) = x ./ sqrt(mean(abs2.(x)))

