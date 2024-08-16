using Test
using DSP
using FiberNlse

const TEST_DIR = dirname(@__FILE__)

for b ∈ ["physical_tests"]
	@info "Testing $(b) bundle"
	for test ∈ readdir(joinpath(dirname(@__FILE__), b), join = true)
		include(test)
	end
end
