using SurrealNumbers
using Documenter

DocMeta.setdocmeta!(SurrealNumbers, :DocTestSetup, :(using SurrealNumbers); recursive=true)

makedocs(;
    modules=[Polylogarithms],
    authors="Matthew Roughan <matthew.roughan@adelaide.edu.au>",
    sitename="SurrealNumbers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://matthew.roughan@adelaide.edu.au.github.io/SurrealNumbers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mroughan/SurrealNumbers.jl",
)
