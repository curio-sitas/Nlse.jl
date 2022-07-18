using Test
using FiberNlse


tests = [
    "soliton",
    "disp_compensation"
]

const testdir = dirname(@__FILE__)


for t in tests
    tp = joinpath(testdir, "$(t).jl")
    include(tp)
end
