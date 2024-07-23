function NL_spm(u, model)
    return 1.0im .* model.γ * model.FFT * (u .* abs2.(u))
end

function NL_spm_self_steepening(u, model)
    η = (1 .+ model.ω / model.ω0)
    return 1.0im .* model.γ * η .* (model.FFT * (u .* abs2.(u)))
end

function NL_spm_self_steepening_raman(u, model)
    η = (1 .+ model.ω / model.ω0)
    IT = abs2.(u)
    RS = model.dt * model.fr * (model.IFFT * ((model.FFT * IT) .* model.RW))
    M = model.FFT * (u .* ((1.0 - model.fr) .* IT .+ RS))
    return 1.0im .* model.γ * η .* M
end

function NL_spm_raman(u, model)
    IT = abs2.(u)
    RS = model.dt * model.fr * (model.IFFT * ((model.FFT * IT) .* model.RW))
    M = model.FFT * (u .* ((1.0 - model.fr) .* IT .+ RS))
    return 1.0im .* model.γ .* M
end

function setup_nonlinearity_term(self_steepening::Bool = true, raman::Bool = true)

    if self_steepening && raman
        return NL_spm_self_steepening_raman
    elseif self_steepening && !raman
        return NL_spm_self_steepening
    elseif !self_steepening && raman
        return NL_spm_raman
    else
        return NL_spm
    end

end