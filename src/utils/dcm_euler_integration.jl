
"""
    dcm_euler_gen(dcm,triple_input)

Performs DCM Euler integration by preparing necessary hemodynamic constants and DCM matrices
and calls the Euler integration function.

# Arguments
- `dcm::DCM`: DCM structure (can be a linear, bilinear on nonlinear DCM)
- `triple_input::Bool`: whether or not to triple the input (to avoid unwanted effects during
convolution with HRF)

# Output
- `y`: BOLD signal
- `x`: neuronal signal
"""
function dcm_euler_gen(dcm::T, triple_input::Bool) where {T<:DCM}
    if !isfile(joinpath(euler_integration_bin, "libdcm_euler_integration.so"))
        @warn "Binary file for Euler integration is missing. Cannot generate synthetic data.
        Please check if gcc is installed and recompile the c code."
        return zeros(dcm.scans, dcm.nr), zeros(dcm.scans, dcm.nr)
    end
    # number of regions
    nr = size(dcm.Ep.A, 1)
    nu = size(dcm.c, 2)
    A = copy(dcm.Ep.A')

    # reserve memory
    B = Array{Float64}(undef, nr, nr, nu)
    D = Array{Float64}(undef, nr, nr, nr)

    # create B and D matrix based on if it's a Linear/Bi-Linear/Non-Linear DCM
    if dcm isa LinearDCM
        B .= zeros(nr, nr, nu)
        D .= zeros(nr, nr, nr)
    elseif dcm isa BiLinearDCM
        B .= permutedims(dcm.Ep.B, (2, 1, 3))
        D .= zeros(nr, nr, nr)
    elseif dcm isa NonLinearDCM
        B .= permutedims(dcm.Ep.B, (2, 1, 3))
        D .= permutedims(dcm.Ep.D, (2, 1, 3))
    end

    # driving inputs
    if triple_input
        u = [dcm.U.u; dcm.U.u; dcm.U.u]
    else
        u = copy(dcm.U.u)
    end
    C = u * (dcm.Ep.C' ./ 16)
    U = u

    # hemodynamic constants
    H = [0.64, 0.32, 2.00, 0.32, 0.32]

    # constants for hemodynamic model
    oxygenExtractionFraction = 0.32 * ones(nr)
    α_inv = 1.0 / H[4]
    τ = H[3] * exp.(dcm.Ep.transit)
    γ = H[2]
    κ = H[1] * exp.(dcm.Ep.decay)
    ϵ = 1.0 * exp.(dcm.Ep.epsilon) * ones(nr)

    # parameter list
    # the last two parameters are for B and D matrix respectively
    # if set to 1 we assume we have bi-linear/non-linear DCM
    if dcm isa LinearDCM
        paramList = vec([dcm.U.dt size(U, 1) nr size(U, 2) 0 0 0])
    elseif dcm isa BiLinearDCM
        paramList = vec([dcm.U.dt size(U, 1) nr size(U, 2) 0 1 0])
    elseif dcm isa NonLinearDCM
        paramList = vec([dcm.U.dt size(U, 1) nr size(U, 2) 0 1 1])
    end
    # neuronal signal and time courses for hemodynamic parameters
    x, _, _, v, q = dcm_euler_integration(
        A, C, U, B, D, oxygenExtractionFraction, α_inv, τ, γ, κ, paramList
    )

    # constants for BOLD signal equation
    relaxationRateSlope = 25.0
    frequencyOffset = 40.3
    oxygenExtractionFraction = 0.4 * ones(nr)
    echoTime = 0.04
    restingVenousVolume = 4.0

    # coefficients of BOLD signal equation
    coefficientK1 = 4.3 * frequencyOffset * echoTime * oxygenExtractionFraction
    coefficientK2 = ϵ .* (relaxationRateSlope * oxygenExtractionFraction * echoTime)
    coefficientK3 = 1.0 .- ϵ'

    # add additional dimension for elementwise multiplication
    coefficientK1 = reshape(coefficientK1, :, 1)
    coefficientK2 = reshape(coefficientK2, :, 1)
    coefficientK3 = reshape(coefficientK3, :, 1)

    Indices = euler_make_indices(dcm)

    y =
        restingVenousVolume * (
            coefficientK1' .* (1.0 .- q[Indices, :]) +
            coefficientK2' .* (1.0 .- (q[Indices, :] ./ v[Indices, :])) +
            coefficientK3' .* (1.0 .- v[Indices, :])
        )

    return y, x
end

"""
    dcm_euler_integration(A,C,U,B,D,rho,alphaInv,tau,gamma,kappa,param)

Wrapper that calls the C code performing the Euler integration.

# Arguments
- `A::Matrix{Float64}`: Endogenous connectivity
- `C::Matrix{Float64}`: Input, represents u*C'
- `U::Matrix{Float64}`: Input matrix
- `B::Array{Float64,3}`: Bi linear connectivity matrix
- `D::Array{Float64,3}`: Non linear connectivity matrix
- `rho::Vector{Float64}`: hemodynamic parameter (one for each region)
- `alphaInv::Float64`: hemodynamic parameter
- `tau::Vector{Float64}`: hemodynamic parameter (one for each region)
- `gamma::Float64`: hemodynamic parameter
- `kappa::Vector{Float64}`: hemodynamic parameter (one for each region)
- `param::Vector{Float64}`: Array of integration specific parameters
                            [Euler step size,
                            total steps,
                            number of regions,
                            number of inputs,
                            subject number]

# Output
- `x::Matrix{Float64}`: neuronal activity
- `s::Matrix{Float64}`: vasodilatory signal
- `f::Matrix{Float64}`: blood flow
- `v::Matrix{Float64}`: blood volume
- `q::Matrix{Float64}`: deoxyhemoglobin content
"""
function dcm_euler_integration(
    A::Matrix{Float64},
    C::Matrix{Float64},
    U::Matrix{Float64},
    B::Array{Float64,3},
    D::Array{Float64,3},
    rho::Vector{Float64},
    alphaInv::Float64,
    tau::Vector{Float64},
    gamma::Float64,
    kappa::Vector{Float64},
    param::Vector{Float64},
)
    nTime = Int(param[2])
    nStates = Int(param[3])

    # allocate memory for the results of the c code
    x_out = Matrix{Cdouble}(undef, nTime, nStates)
    s_out = Matrix{Cdouble}(undef, nTime, nStates)
    f1_out = Matrix{Cdouble}(undef, nTime, nStates)
    v1_out = Matrix{Cdouble}(undef, nTime, nStates)
    q1_out = Matrix{Cdouble}(undef, nTime, nStates)

    @ccall abspath(
        joinpath(euler_integration_bin, "libdcm_euler_integration.so")
    ).dcm_euler_integration(
        A::Ref{Cdouble},
        C::Ref{Cdouble},
        U::Ref{Cdouble},
        B::Ref{Cdouble},
        D::Ref{Cdouble},
        rho::Ref{Cdouble},
        alphaInv::Cdouble,
        tau::Ref{Cdouble},
        gamma::Cdouble,
        kappa::Ref{Cdouble},
        param::Ref{Cdouble},
        x_out::Ref{Cdouble},
        s_out::Ref{Cdouble},
        f1_out::Ref{Cdouble},
        v1_out::Ref{Cdouble},
        q1_out::Ref{Cdouble},
    )::Cvoid

    return x_out, s_out, f1_out, v1_out, q1_out
end

function euler_make_indices(dcm::T) where {T<:DCM}
    # length of input
    L = size(dcm.U.u, 1)

    # In rDCM code (matlab) a vector of ones is created but only the first value of delay is used.
    delay = 1 # TODO: make this adjustable by user, specify in DCM or as option

    # create index array
    Indices = collect(1:L)

    # get the indices that coincide wiht the data timepoints
    idx = Vector{Bool}(undef, L)
    idx .= false
    idx[delay:Int(floor(dcm.Y.dt / dcm.U.dt)):end] .= true

    # asign those timings
    return Indices[idx]
end
