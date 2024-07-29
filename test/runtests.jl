using Aqua: Aqua
using JET
using JuliaFormatter: JuliaFormatter
using MAT: matopen
using Test
using RegressionDynamicCausalModeling
# using Scratch
# using LazyArtifacts
using Random
using SparseArrays: spzeros

@testset verbose = true "rDCM tests" begin
    @testset "Code formatting" begin
        @test JuliaFormatter.format(
            RegressionDynamicCausalModeling; verbose=false, overwrite=false
        )
    end

    @testset "Code linting" begin
        JET.test_package(RegressionDynamicCausalModeling; target_defined_modules=true)
    end

    @testset "Code quality" begin
        Aqua.test_all(
            RegressionDynamicCausalModeling;
            ambiguities=false,
            deps_compat=(check_extras=false,),
        )
    end

    for name in (
        "utils",
        "constructors",
        "simulate_dcm",
        "prior",
        "regressors",
        "model_inversion",
        "predict",
    )
        include("test_$name.jl")
    end
end
