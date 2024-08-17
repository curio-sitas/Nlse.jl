function raman_linagrawaal(fr = 0.18, τl = 32e-15,
	τvib = 12.2e-15)

	function time_response(t)
		hr = 0.0 .+ (t .>= 0) .* ((τvib^2 + τl^2) / τvib / τl^2 * exp.(-t / τl) .* sin.(t / τvib))
		return hr
	end

	return RamanModel(fr, time_response)

end
