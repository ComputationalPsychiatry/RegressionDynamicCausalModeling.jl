function get_precision_component(N::Int, TR::Float64, aᵣᵣ::Float64)
    diag0 = zeros(N) .+ 2/TR^2 .+ aᵣᵣ^2
    diag1 = zeros(N-1) .- 1/TR^2 .+ aᵣᵣ^2

    return SymTridiagonal(diag0, diag1) \ Matrix{Float64}(I, N, N)
end
