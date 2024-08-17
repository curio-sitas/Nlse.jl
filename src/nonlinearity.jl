
function NL_spm(u, model)
	return 1.0im .* model.γ * model.fftp * (u .* abs2.(u))
end

function NL_spm_self_steepening(u, model)
	η = (1 .+ model.ω / model.ω0)
	return 1.0im .* model.γ * η .* (model.fftp * (u .* abs2.(u)))
end

function NL_spm_self_steepening_raman(u, model)
	η = (1 .+ model.ω / model.ω0)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	M = model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS))
	return 1.0im .* model.γ * η .* M
end

function NL_spm_raman(u, model)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	M = model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS))
	return 1.0im .* model.γ .* M
end