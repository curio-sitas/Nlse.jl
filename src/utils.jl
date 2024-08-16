

function combine(sol1::Solution, sol2::Solution)

	t = sol1.t
	f = sol1.f
	At = vcat(sol1.At, sol2.At)
	Af = vcat(sol1.Af, sol2.Af)
	z = vcat(sol1.z, sol2.z .+ sol1.z[end])
	return Solution(z, t, f, At, Af)

end

function combine(sols::Vector{Solution})
	sol = sols[1]
	for soli ∈ sols[2:end]
		sol = combine(sol, soli)
	end
	sol
end

