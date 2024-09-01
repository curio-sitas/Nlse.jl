using FiberNlse
using Documenter

makedocs(;
	sitename = "FiberNlse.jl",
	modules = [FiberNlse],
	authors = "curio-sitas <brian.sinquin@gmail.com> and contributors",
	format = Documenter.HTML(;
		canonical = "https://brian-sinquin.github.io/FiberNlse.jl",
		edit_link = "main",
		assets = String[],
	),
	pages = [
		"Home" => "index.md",
		#joinpath.("examples", filter(x -> endswith(x, ".md"), readdir(MD_OUTPUT))),
		"API" => "api.md",
	],
)

deploydocs(;
	repo = "github.com/brian-sinquin/FiberNlse.jl",
	push_preview = true,
)
