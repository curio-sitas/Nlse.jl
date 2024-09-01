using FiberNlse
using Documenter

DocMeta.setdocmeta!(FiberNlse, :DocTestSetup, :(using FiberNlse); recursive = true)

makedocs(;
	modules = [FiberNlse],
	authors = "curio-sitas <brian.sinquin@gmail.com> and contributors",
	sitename = "FiberNlse.jl",
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
