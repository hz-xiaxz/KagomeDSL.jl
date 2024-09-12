using KagomeDSL
using Documenter

DocMeta.setdocmeta!(KagomeDSL, :DocTestSetup, :(using KagomeDSL); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
  file for
  file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [KagomeDSL],
    authors = "hzxiaxz <hzxiaxz@gmail.com>",
    repo = "https://github.com/hz-xiaxz/KagomeDSL.jl/blob/{commit}{path}#{line}",
    sitename = "KagomeDSL.jl",
    format = Documenter.HTML(; canonical = "https://hz-xiaxz.github.io/KagomeDSL.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/hz-xiaxz/KagomeDSL.jl")
