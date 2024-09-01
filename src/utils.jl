"""
Combines two different [Solution](@ref) structs to an unique one.
"""
function combine(sol1::Solution, sol2::Solution)

	t = sol1.t
	f = sol1.f
	At = vcat(sol1.At, sol2.At)
	Af = vcat(sol1.Af, sol2.Af)
	z = vcat(sol1.z, sol2.z .+ sol1.z[end])
	Solution(z, t, f, At, Af)

end


"""
Returns the time and frequency domain output signal of a [Solution](@ref)
"""
function output(sol::Solution)
	return sol.At[end, :], sol.Af[end, :]
end

"""
Returns the time and frequency domain input signal of a [Solution](@ref)
"""
function input(sol::Solution)
	return sol.At[1, :], sol.Af[1, :]
end

"""
Combines a set of [Solution](@ref) structs to an unique one.
"""
function combine(sols::Vector{Solution})
	sol = sols[1]
	for soli ∈ sols[2:end]
		sol = combine(sol, soli)
	end
	sol
end

function chirp(x)
	@error "Not implemented!"
end
