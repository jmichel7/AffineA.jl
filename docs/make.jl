using Documenter, AffineA, Chevie

DocMeta.setdocmeta!(AffineA, :DocTestSetup, :(using AffineA); recursive=true)

makedocs(;
    modules=[AffineA],
    authors="Jean Michel <jean.michel@imj-prg.fr> and contributors",
    sitename="AffineA.jl",
    format=Documenter.HTML(;
        canonical="https://jmichel7.github.io/AffineA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=:missing_docs,
)

deploydocs(;
    repo="github.com/jmichel7/AffineA.jl",
    devbranch="main",
)
