
#TODO!: inplace functions

function _spm(u, model)
	return 1.0im .* model.γ * (model.fftp * (u .* abs2.(u)))
end

function _spm_self_steepening(u, model)
	return 1.0im .* model.γ * (1 .+ model.ω / model.ω0) .* (model.fftp * (u .* abs2.(u)))
end

function _spm_self_steepening_raman(u, model)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	return 1.0im .* model.γ * (1 .+ model.ω / model.ω0) .* (model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS)))
end

function _spm_raman(u, model)
	IT = abs2.(u)
	RS = model.dt * model.fr * (model.ifftp * ((model.fftp * IT) .* model.raman_freq_response))
	return 1.0im .* model.γ .* (model.fftp * (u .* ((1.0 - model.fr) .* IT .+ RS)))
end

function choose_nonlinear_term(self_steepening::Bool = true, raman::Bool = true)

	if self_steepening && raman
		return _spm_self_steepening_raman
	elseif self_steepening && !raman
		return _spm_self_steepening
	elseif !self_steepening && raman
		return _spm_raman
	else
		return _spm
	end

end
