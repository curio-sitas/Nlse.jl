using FiberNlse
using DSP

using DataFrames
using CSV
using Test


for bundle ∈ ["api_tests", "physical_tests", "comparison_tests"]
	@info "Testing $(bundle) bundle"
	test_files = filter(!isdir, readdir(joinpath(dirname(@__FILE__), bundle), join = true))
	for test ∈ test_files
		include(test)
	end
end
