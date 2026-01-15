
"""
    dcm_euler_gen(dcm)

Performs DCM Euler integration by preparing necessary hemodynamic constants and DCM matrices
and calls the Euler integration function.

# Arguments
- `dcm::DCM`: DCM structure (can be a linear, bilinear on nonlinear DCM)

# Output
- `y`: BOLD signal
- `x`: neuronal signal
"""
function dcm_euler_gen(dcm::T) where {T<:DCM}

    # number of regions
    nr = size(dcm.a, 1)
    A = copy(dcm.Ep.A')

    # create B and D matrix based on if it's a Linear/Bi-Linear/Non-Linear DCM
    B, D = get_matrices(dcm)

    # driving inputs
    C = dcm.U.u * (dcm.Ep.C' ./ 16)
    U = copy(dcm.U.u)

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
    paramList = get_param_list(dcm, U, nr)

    # neuronal signal and time courses for hemodynamic parameters
    if isfile(joinpath(euler_integration_bin, "libdcm_euler_integration.so"))
        x, _, _, v, q = dcm_euler_integration_c(
            A, C, U, B, D, oxygenExtractionFraction, α_inv, τ, γ, κ, paramList
        )
    else
        x, _, _, v, q = dcm_euler_integration_jl(
            A,
            C,
            U,
            B,
            D,
            oxygenExtractionFraction,
            α_inv,
            τ,
            γ,
            κ,
            dcm.U.dt,
            size(U, 1),
            nr,
            size(U, 2),
            Bool(paramList[6]),
            Bool(paramList[7]),
        )
    end

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

function get_matrices(dcm::LinearDCM)
    nr = size(dcm.a, 1)
    nu = size(dcm.c, 2)
    B = zeros(nr, nr, nu)
    D = zeros(nr, nr, nr)

    return B, D
end

function get_matrices(dcm::BiLinearDCM)
    nr = size(dcm.a, 1)
    B = permutedims(dcm.Ep.B, (2, 1, 3))
    D = zeros(nr, nr, nr)

    return B, D
end

function get_matrices(dcm::NonLinearDCM)
    B = permutedims(dcm.Ep.B, (2, 1, 3))
    D = permutedims(dcm.Ep.D, (2, 1, 3))

    return B, D
end

function get_param_list(dcm::LinearDCM, U, nr)
    return vec([dcm.U.dt size(U, 1) nr size(U, 2) 0 0 0])
end

function get_param_list(dcm::BiLinearDCM, U, nr)
    return vec([dcm.U.dt size(U, 1) nr size(U, 2) 0 1 0])
end

function get_param_list(dcm::NonLinearDCM, U, nr)
    return vec([dcm.U.dt size(U, 1) nr size(U, 2) 0 1 1])
end

"""
    dcm_euler_integration_c(A,C,U,B,D,rho,alphaInv,tau,gamma,kappa,param)

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
function dcm_euler_integration_c(
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

    @ccall abspath(joinpath(euler_integration_bin, "libdcm_euler_integration.so")).dcm_euler_integration(
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

"""
    dcm_euler_integration_jl(A,C,U,B,D,rho,alphaInv,tau,gamma,kappa,param)

Perform Euler integration of DCM.

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
- `timeStep::Float64`: Euler step size
- `nTime::Int`: total steps
- `nStates::Int`: number of regions
- `nInputs::Int`: number of inputs
- `dcmTypeB::Bool`: is it a Bilinear DCM
- `dcmTypeD::Bool`: is it a non-linear DCM

# Output
- `x::Matrix{Float64}`: neuronal activity
- `s::Matrix{Float64}`: vasodilatory signal
- `f::Matrix{Float64}`: blood flow
- `v::Matrix{Float64}`: blood volume
- `q::Matrix{Float64}`: deoxyhemoglobin content
"""
function dcm_euler_integration_jl(
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
    timeStep::Float64,
    nTime::Int,
    nStates::Int,
    nInputs::Int,
    dcmTypeB::Bool,
    dcmTypeD::Bool,
)

    # Initialize the dynamical system to resting state values
    x_out = zeros(nTime, nStates)
    s_out = zeros(nTime, nStates)
    f1_out = ones(nTime, nStates)
    v1_out = ones(nTime, nStates)
    q1_out = ones(nTime, nStates)
    old_f1 = zeros(nStates)
    old_v1 = zeros(nStates)
    old_q1 = zeros(nStates)

    # Euler's integration steps
    for iStep in 1:(nTime - 1)
        # for each region
        for jState in 1:nStates
            # update x (neuronal signal)
            x_out[iStep + 1, jState] =
                x_out[iStep, jState] +
                timeStep * (A[:, jState]' * x_out[iStep, :] + C[iStep, jState])

            if dcmTypeB
                # B matrix update
                for kIter in 1:nInputs
                    x_out[iStep + 1, jState] +=
                        timeStep *
                        (U[iStep, kIter] * x_out[iStep, :]' * B[:, jState, kIter])
                end
            end

            if dcmTypeD
                # D matrix update
                for kIter in 1:nStates
                    x_out[iStep + 1, jState] +=
                        timeStep *
                        (x_out[iStep, kIter] * (x_out[iStep, :]' * D[:, jState, kIter]))
                end
            end

            # update s (vasodilatory signal)
            s_out[iStep + 1, jState] =
                s_out[iStep, jState] +
                timeStep * (
                    x_out[iStep, jState] - kappa[jState] * s_out[iStep, jState] -
                    gamma * (f1_out[iStep, jState] - 1.0)
                )

            # update f (blood flow)
            old_f1[jState] += timeStep * (s_out[iStep, jState] / f1_out[iStep, jState])

            # update v (blood volume)
            tmp1 = v1_out[iStep, jState]^alphaInv
            tmp2 = (1.0 - (1.0 - rho[jState])^(1.0 / f1_out[iStep, jState])) / rho[jState]

            old_v1[jState] +=
                timeStep *
                ((f1_out[iStep, jState] - tmp1) / (tau[jState] * v1_out[iStep, jState]))

            # update q (deoxyhemoglobin content)
            old_q1[jState] +=
                timeStep * (
                    (
                        f1_out[iStep, jState] * tmp2 -
                        tmp1 * q1_out[iStep, jState] / v1_out[iStep, jState]
                    ) / (tau[jState] * q1_out[iStep, jState])
                )

            # tracking the exponential values
            f1_out[iStep + 1, jState] = exp(old_f1[jState])
            v1_out[iStep + 1, jState] = exp(old_v1[jState])
            q1_out[iStep + 1, jState] = exp(old_q1[jState])
        end
    end

    return x_out, s_out, f1_out, v1_out, q1_out
end
