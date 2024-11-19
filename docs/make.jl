using Documenter
using RegressionDynamicCausalModeling
using Literate

open(joinpath(joinpath(@__DIR__, "src"), "index.md"), "w") do io
    println(
        io,
        """
        ```@meta
        EditURL = "https://github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl/blob/main/README.md"
        ```
        """,
    )
    for line in eachline(joinpath(dirname(@__DIR__), "README.md"))
        println(io, line)
    end
end

examples_jl_path = joinpath(dirname(@__DIR__), "examples")
examples_md_path = joinpath(@__DIR__, "src", "examples")

for file in readdir(examples_md_path)
    if endswith(file, ".md")
        rm(joinpath(examples_md_path, file))
    end
end

for file in readdir(examples_jl_path)
    Literate.markdown(
        joinpath(examples_jl_path, file),
        examples_md_path;
        documenter=true,
        #codefence =  "```text" => "```"
    )
end

pages = [
    "Home" => "index.md",
    "API reference" => "api.md",
    "Tutorials" => [
        "SPM compatibility" => joinpath("examples", "SPM_compat.md"),
        "Construct DCM" => joinpath("examples", "Specify_DCM.md"),
        "Parameter estimation" => joinpath("examples", "Invert_rDCM.md"),
    ],
]

fmt = Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    repolink="https://github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl",
    canonical="https://ComputationalPsychiatry.github.io/RegressionDynamicCausalModeling.jl",
    assets=String[],
    collapselevel=1,
    description="A Julia package for estimating effective (i.e., directed) connectivity in large (whole-brain) networks",
    size_threshold_ignore=[joinpath("examples", "Invert_rDCM.md")],
    #sidebar_sitename=false,
)

makedocs(;
    sitename="RegressionDynamicCausalModeling",
    modules=[RegressionDynamicCausalModeling],
    authors="Imre Kertesz",
    format=fmt,
    pages=pages,
    checkdocs=:exports,
    pagesonly=true,
)

deploydocs(;
    repo="github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl.git",
    devbranch="main",
)
