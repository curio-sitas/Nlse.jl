using Documenter, Literate
using FiberNlse

EXAMPLES = joinpath("../", @__DIR__, "examples")
MD_OUTPUT = joinpath(@__DIR__, "src", "examples")
ispath(MD_OUTPUT) && rm(MD_OUTPUT; recursive = true)


for file in readdir(EXAMPLES; join = true)
	endswith(file, ".jl") || continue
	Literate.markdown(file, MD_OUTPUT; documenter = true)
end

makedocs(
	sitename = "FiberNlse",
	format = Documenter.HTML(),
	modules = [FiberNlse],
	authors = "Brian Sinquin",
	pages = [
		"Home" => "index.md",
		"Theoretical Background" => "theory.md",
		"User Guide" => "userguide.md",
		"Examples" =>
			joinpath.("examples", filter(x -> endswith(x, ".md"), readdir(MD_OUTPUT))),
		"API" => "api.md",
	],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(;
	repo = "github.com/brian-sinquin/FiberNlse.jl.git",
	target = "build",
	push_preview = true,
)
