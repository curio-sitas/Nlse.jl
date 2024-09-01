using Test

using FiberNlse
using DSP

using DataFrames
using CSV


for bundle ∈ ["comparison_tests"]
	@info "Testing $(bundle) bundle"
	for test ∈ readdir(joinpath(dirname(@__FILE__), bundle), join = true)
		include(test)
	end
end
