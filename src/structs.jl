abstract type AbstractInvParam end
abstract type ModelOutput end

"""
    DCM

Abstract supertype for a DCM.
"""
abstract type DCM end

"""
    RDCM

Abstract supertype for a rigid or sparse rDCM.
"""
abstract type RDCM end
abstract type TrueParam end

#---------------------------------------------------------------------------
# definition of structs
#---------------------------------------------------------------------------
"""
$(TYPEDEF)

Settings for model inversion specific to rigid rDCM (fixed network architecture).

# Fields

$(TYPEDFIELDS)
"""
struct RigidInversionParams <: AbstractInvParam
    "Maximum number of iterations per region"
    maxIter::Int
    "Tolerance for convergence"
    tol::Float64

    # inner constructor with sanity checks
    function RigidInversionParams(maxIter, tol)
        if maxIter ≤ 0 || tol ≤ 0.0
            error("Invalid parameters.")
        end
        return new(maxIter, tol)
    end
end

"""
$(TYPEDEF)

Settings for model inversion specific to sparse rDCM.

# Fields

$(TYPEDFIELDS)
"""
struct SparseInversionParams <: AbstractInvParam
    "Maximum number of iterations per region"
    maxIter::Int
    "Tolerance for convergence"
    tol::Float64
    "Number of reruns"
    reruns::Int
    "Whether or not to estimate sparsity for C matrix. If true, the Bernoulli posteriors for the C matrix are not pruned."
    restrictInputs::Bool

    # inner constructor with sanity checks
    function SparseInversionParams(maxIter, tol, reruns, restrictInputs)
        if maxIter ≤ 0 || tol ≤ 0.0 || reruns ≤ 0
            error("Invalid inversion parameters.")
        end
        return new(maxIter, tol, reruns, restrictInputs)
    end
end

Base.@kwdef mutable struct RegressorParams
    padding::Bool = false
    #filterStr::Int = 0
end

"""
$(TYPEDEF)

Settings for model inversion.

# Fields

$(TYPEDFIELDS)
"""
struct Options{T<:AbstractInvParam}
    "Model specific inversion settings."
    invParams::T
    "Verbosity during inversion (0,1 or 2)"
    verbose::Int
    "Whether or not synthetic data is used."
    synthetic::Bool
    "Used for testing"
    testing::Bool
    "Random number generator"
    rng::AbstractRNG

    # inner constructor with sanity checks
    function Options(invParams, verbose, synthetic, testing, rng)
        if verbose < 0
            error("Verbosity level needs to be an integer between 0 and 2.")
        end
        return new{typeof(invParams)}(invParams, verbose, synthetic, testing, rng)
    end
end

"""
$(TYPEDEF)

Output after inversion of [`RigidRdcm`](@ref RigidRdcm).

# Fields

$(TYPEDFIELDS)
"""
struct RigidOutput <: ModelOutput
    "Negative free energy of the whole model"
    F::Float64
    "Region specific negative free energy"
    F_r::Vector{Float64}
    "Number of iterations per region until convergence"
    iter_all::Vector{Int}
    "Posterior shape parameter for noise"
    a_all::Vector{Float64}
    "Posterior rate parameter for noise"
    b_all::Vector{Float64}
    "Posterior mean for connectivity parameters"
    m_all::Matrix{Float64}
    "Posterior covariance for connectivity parameters"
    Σ_all::Vector{SparseMatrixCSC{Float64,Int}}
    "Inversion method"
    inversion::String

    # inner constructor with sanity checks
    function RigidOutput(F, F_r, iter_all, a_all, b_all, m_all, Σ_all, inversion)
        if sum(F_r) ≠ F
            error(
                "Sum of region-wise neg. free energies don't sum up to overall neg. free energy.",
            )
        end
        if any(a_all .≤ 0.0) || any(b_all .≤ 0.0)
            error("Found invalid values of the posterior Gamma distribution.")
        end
        if any(iter_all .≤ 0)
            error("Invalid number of iterations")
        end
        if length(a_all) ≠ length(b_all) ||
            length(b_all) ≠ size(m_all, 1) ||
            size(m_all, 1) ≠ size(Σ_all, 1)
            error("Inconsistent number of regions.")
        end
        nr = length(F_r)
        idx = m_all .≠ 0.0
        for r in 1:nr
            idxᵣ = idx[r, :]
            if !isposdef(Σ_all[r][idxᵣ, idxᵣ])
                error("One or more covariance matrices are not positive definite.")
            end
        end

        return new(F, F_r, iter_all, a_all, b_all, m_all, Σ_all, inversion)
    end
end

"""
$(TYPEDEF)

Output after inversion of [`SparseRdcm`](@ref SparseRdcm).

# Fields

$(TYPEDFIELDS)
"""
struct SparseOutput <: ModelOutput
    "Negative free energy of the whole model"
    F::Float64
    "Region specific negative free energy"
    F_r::Vector{Float64}
    "Number of iterations per region until convergence"
    iter_all::Vector{Int} # number of iterations per region until convergence
    "Posterior shape parameter for noise"
    a_all::Vector{Float64}
    "Posterior rate parameter for noise"
    b_all::Vector{Float64}
    "Posterior mean for connectivity parameters"
    m_all::Matrix{Float64}
    "Posterior covariance for connectivity parameters"
    Σ_all::Vector{SparseMatrixCSC{Float64,Int}}
    "Posterior for binary indicator variables"
    z_all::Matrix{Float64}
    "Inversion method"
    inversion::String

    # inner constructor with sanity checks
    function SparseOutput(F, F_r, iter_all, a_all, b_all, m_all, Σ_all, z_all, inversion)
        if sum(F_r) ≠ F
            error(
                "Sum of region-wise neg. free energies don't sum up to overall neg. free energy.",
            )
        end
        if any(a_all .≤ 0.0) || any(b_all .≤ 0.0)
            error("Found invalid values of the posterior Gamma distribution.")
        end
        if any(iter_all .≤ 0)
            error("Invalid number of iterations")
        end
        if any(z_all .< 0.0) || any(z_all .> 1.0)
            error("Invalid probabilities in posterior Bernoulli.")
        end
        if length(a_all) ≠ length(b_all) ||
            length(b_all) ≠ size(m_all, 1) ||
            size(m_all, 1) ≠ size(Σ_all, 1) ||
            size(Σ_all, 1) ≠ size(z_all, 1)
            error("Inconsistent number of regions.")
        end
        nr = length(F_r)
        idx = m_all .≠ 0.0
        for r in 1:nr
            idxᵣ = idx[r, :]
            if !isposdef(Σ_all[r][idxᵣ, idxᵣ])
                error("One or more covariance matrices are not positive definite.")
            end
        end

        return new(F, F_r, iter_all, a_all, b_all, m_all, Σ_all, z_all, inversion)
    end
end

"""
$(TYPEDEF)

Input structure for a DCM.

# Fields

$(TYPEDFIELDS)
"""
struct InputU
    "Driving input"
    u::Matrix{Float64}
    "Sampling interval"
    dt::Float64
    "Input names"
    name::Vector{String}

    # inner constructor with sanity checks
    function InputU(u, dt, name)
        if size(u, 2) ≠ length(name)
            error("Size of u and name vector don't match.")
        end
        if dt ≤ 0.0
            error("Sampling rate must be positive.")
        end
        return new(u, dt, name)
    end
end

"""
$(TYPEDEF)

Information about confounds.

# Fields

$(TYPEDFIELDS)
"""
struct Confound
    "Confound matrix of size number of datapoints times number of confounds."
    X0::Union{Vector{Float64},Matrix{Float64}} # confounds
    "Confound names"
    name::Vector{String}                       # confound names

    # inner constructor with sanity checks
    function Confound(X0, name)
        if size(X0, 2) ≠ length(name)
            error("Size of X0 and name vector don't match.")
        end
        return new(X0, name)
    end
end

"""
$(TYPEDEF)

Struct for BOLD signal.

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct BoldY
    "BOLD signal"
    y::Union{Matrix{Float64},Nothing}
    "Sampling interval"
    dt::Float64
    "Brain region names"
    name::Union{Vector{String},Nothing}

    # inner constructor with sanity checks
    function BoldY(y, dt, name)
        if !isnothing(y) && !isnothing(name) && size(y, 2) ≠ length(name)
            error("Size of y and name vector don't match.")
        end
        if dt ≤ 0.0
            error("Sampling rate must be positive.")
        end
        return new(y, dt, name)
    end
end

struct TrueParamLinear <: TrueParam
    A::Matrix{Float64}
    C::Matrix{Float64}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64

    # inner constructor with sanity checks
    function TrueParamLinear(A, C, transit, decay, epsilon)
        if size(A, 1) ≠ size(C, 1)
            error("Size of A and C matrix don't match.")
        end
        if size(transit, 1) ≠ size(decay, 1) || size(transit, 1) ≠ size(A, 1)
            error("Inconsistent number of regions.")
        end
        return new(A, C, transit, decay, epsilon)
    end
end

struct TrueParamBiLinear <: TrueParam
    A::Matrix{Float64}
    B::Array{Float64,3}
    C::Matrix{Float64}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64

    # inner constructor with sanity checks
    function TrueParamBiLinear(A, B, C, transit, decay, epsilon)
        if size(A, 1) ≠ size(C, 1) || size(B, 1) ≠ size(A, 1) || size(B, 3) ≠ size(C, 2)
            error("Size of A, B or C matrix don't match.")
        end
        if size(transit, 1) ≠ size(decay, 1) || size(transit, 1) ≠ size(A, 1)
            error("Inconsistent number of regions.")
        end
        return new(A, B, C, transit, decay, epsilon)
    end
end

struct TrueParamNonLinear <: TrueParam
    A::Matrix{Float64}
    B::Array{Float64,3}
    C::Matrix{Float64}
    D::Array{Float64,3}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64

    # inner constructor with sanity checks
    function TrueParamNonLinear(A, B, C, D, transit, decay, epsilon)
        if size(A, 1) ≠ size(C, 1) ||
            size(B, 1) ≠ size(A, 1) ||
            size(B, 3) ≠ size(C, 2) ||
            size(D, 1) ≠ size(A, 1) ||
            size(D, 1) ≠ size(D, 2) ||
            size(D, 2) ≠ size(D, 3)
            error("Size of A, B or C matrix don't match.")
        end
        if size(transit, 1) ≠ size(decay, 1) || size(transit, 1) ≠ size(A, 1)
            error("Inconsistent number of regions.")
        end
        return new(A, B, C, D, transit, decay, epsilon)
    end
end

"""
$(TYPEDEF)

Representation of a linear DCM.

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct LinearDCM <: DCM
    "Binary indicator matrix for endogenous connectivity"
    a::BitMatrix
    "Binary indicator matrix for driving inputs"
    c::BitMatrix
    "number of data points per region"
    scans::Int
    "number of regions in network"
    const nr::Int
    "input structure with information about driving input"
    U::Union{InputU,Nothing}
    "data structure containing BOLD signal"
    Y::Union{BoldY,Nothing}
    "connectivity parameters (A and C matrix)"
    Ep::TrueParamLinear
    "confound structure"
    Conf::Union{Confound,Nothing}

    @doc """
    $(SIGNATURES)

    Constructor with sanity checks for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
    """
    function LinearDCM(a, c, scans, nr, U, Y, Ep, Conf)
        if nr ≠ size(a, 1) ||
            size(c, 1) ≠ nr ||
            (!isnothing(Y) && !isnothing(Y.y) && size(Matrix{Float64}(Y.y), 2) ≠ nr) ||
            size(Ep.A, 1) ≠ nr
            error("Inconsistent number of regions.")
        end
        if (!isnothing(Conf) && !isnothing(U) && size(Conf.X0, 1) ≠ size(U.u, 1))
            error("Confound matrix size and input matrix size don't match.")
        end
        if scans ≤ 0
            error("Invalid number of scans.")
        end
        if !isnothing(Y) && !isnothing(U)
            y = Y.y # for JET not to throw error
            if !isnothing(y)
                r_dt = 1
                try
                    r_dt = Int(Y.dt / U.dt)
                catch
                    error(
                        "The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                  of the input U (u_dt). Cannot proceed.",
                    )
                end
                if size(y, 1) ≠ size(U.u, 1) / r_dt
                    error("Length of BOLD signal and driving input u is inconsistent.")
                end
            end
        end

        return new(a, c, scans, nr, U, Y, Ep, Conf)
    end
end

"""
$(TYPEDEF)

Representation of a bi-linear DCM.

# Fields

$(TYPEDFIELDS)

!!! warning
    While this package allows to simulate data for bi-linear DCMs, the current rDCM implementation can only invert linear DCMs.

"""
Base.@kwdef mutable struct BiLinearDCM <: DCM
    "Binary indicator matrix for endogenous connectivity"
    a::BitMatrix
    "Binary indicator matrix for bi-linear dynamics."
    b::BitArray{3}
    "Binary indicator matrix for driving inputs"
    c::BitMatrix
    "number of data points per region"
    scans::Int
    "number of regions in network"
    const nr::Int
    "input structure with information about driving input"
    U::Union{InputU,Nothing}
    "data structure containing BOLD signal"
    Y::Union{BoldY,Nothing}
    "connectivity parameters (A, B and C matrix)"
    Ep::TrueParamBiLinear
    "confound structure"
    Conf::Union{Confound,Nothing}

    @doc """
    $(SIGNATURES)

    Constructor with sanity checks for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
    """
    function BiLinearDCM(a, b, c, scans, nr, U, Y, Ep, Conf)
        if nr ≠ size(a, 1) ||
            size(c, 1) ≠ nr ||
            (!isnothing(Y) && !isnothing(Y.y) && size(Y.y, 2) ≠ nr) ||
            size(Ep.A, 1) ≠ nr ||
            size(Ep.B, 1) ≠ nr ||
            size(b, 1) ≠ nr
            error("Inconsistent number of regions.")
        end
        if size(c, 2) ≠ size(b, 3)
            error("Number of inputs don't match.")
        end
        if (!isnothing(Conf) && !isnothing(U) && size(Conf.X0, 1) ≠ size(U.u, 1))
            error("Confound matrix size and input matrix size don't match.")
        end
        if scans ≤ 0
            error("Invalid number of scans.")
        end
        if !isnothing(Y) && !isnothing(Y.y) && !isnothing(U)
            r_dt = 1
            try
                r_dt = Int(Y.dt / U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
            if size(Y.y, 1) ≠ size(U.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsistent.")
            end
        end

        return new(a, b, c, scans, nr, U, Y, Ep, Conf)
    end
end

"""
$(TYPEDEF)

Representation of a non-linear DCM.

# Fields

$(TYPEDFIELDS)

!!! warning
    While this package allows to simulate data for non-linear DCMs, the current rDCM implementation can only invert linear DCMs.

"""
Base.@kwdef mutable struct NonLinearDCM <: DCM
    "Binary indicator matrix for endogenous connectivity"
    a::BitMatrix
    "Binary indicator matrix for bi-linear dynamics."
    b::BitArray{3}
    "Binary indicator matrix for driving inputs"
    c::BitMatrix
    "Binary indicator matrix for non-linear dynamics"
    d::BitArray{3}
    "number of data points per region"
    scans::Int
    "number of regions in network"
    const nr::Int
    "input structure with information about driving input"
    U::Union{InputU,Nothing}
    "data structure containing BOLD signal"
    Y::Union{BoldY,Nothing}
    "connectivity parameters (A, B, C and D matrix)"
    Ep::TrueParamNonLinear
    "confound structure"
    Conf::Union{Confound,Nothing}

    @doc """
    $(SIGNATURES)

    Constructor with sanity checks for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
    """
    function NonLinearDCM(a, b, c, d, scans, nr, U, Y, Ep, Conf)
        if nr ≠ size(a, 1) ||
            size(c, 1) ≠ nr ||
            (!isnothing(Y) && !isnothing(Y.y) && size(Y.y, 2) ≠ nr) ||
            size(Ep.A, 1) ≠ nr ||
            size(Ep.B, 1) ≠ nr ||
            size(b, 1) ≠ nr ||
            size(d, 1) ≠ nr ||
            size(d, 1) ≠ size(d, 2) ||
            size(d, 2) ≠ size(d, 3) ||
            size(d, 1) ≠ size(Ep.D, 1)
            error("Inconsistent number of regions.")
        end
        if size(c, 2) ≠ size(b, 3)
            error("Number of inputs don't match.")
        end
        if (!isnothing(Conf) && !isnothing(U) && size(Conf.X0, 1) ≠ size(U.u, 1))
            error("Confound matrix size and input matrix size don't match.")
        end
        if scans ≤ 0
            error("Invalid number of scans.")
        end
        if !isnothing(Y) && !isnothing(Y.y) && !isnothing(U)
            r_dt = 1
            try
                r_dt = Int(Y.dt / U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
            if size(Y.y, 1) ≠ size(U.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsistent.")
            end
        end

        return new(a, b, c, d, scans, nr, U, Y, Ep, Conf)
    end
end

#Base.@kwdef mutable struct priorMean
struct priorMeanLinear
    A::Matrix{Float64}
    C::Matrix{Float64}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64

    # inner constructor with sanity checks
    function priorMeanLinear(A, C, transit, decay, epsilon)
        if size(A, 1) ≠ size(C, 1) ||
            size(transit, 1) ≠ size(decay, 1) ||
            size(transit, 1) ≠ size(A, 1)
            error("Inconsistent number of regions.")
        end
        return new(A, C, transit, decay, epsilon)
    end
end

struct priorMeanBiLinear
    A::Matrix{Float64}
    B::Array{Float64,3}
    C::Matrix{Float64}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64

    # inner constructor with sanity checks
    function priorMeanBiLinear(A, B, C, transit, decay, epsilon)
        if size(A, 1) ≠ size(B, 1) ||
            size(B, 2) ≠ size(C, 1) ||
            size(B, 1) ≠ size(B, 2) ||
            size(transit, 1) ≠ size(decay, 1) ||
            size(transit, 1) ≠ size(A, 1)
            error("Inconsistent number of regions.")
        end
        return new(A, B, C, transit, decay, epsilon)
    end
end

struct priorMeanNonLinear
    A::Matrix{Float64}
    B::Array{Float64,3}
    C::Matrix{Float64}
    D::Array{Float64,3}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64

    # inner constructor with sanity checks
    function priorMeanNonLinear(A, B, C, D, transit, decay, epsilon)
        if size(A, 1) ≠ size(B, 1) ||
            size(B, 2) ≠ size(C, 1) ||
            size(B, 1) ≠ size(B, 2) ||
            size(D, 1) ≠ size(A, 1) ||
            size(D, 1) ≠ size(D, 2) ||
            size(D, 2) ≠ size(D, 3) ||
            size(transit, 1) ≠ size(decay, 1) ||
            size(transit, 1) ≠ size(A, 1)
            error("Inconsistent number of regions.")
        end
        return new(A, B, C, D, transit, decay, epsilon)
    end
end

#Base.@kwdef mutable struct priorCov
struct priorCovLinear
    A::Matrix{Float64}
    C::Matrix{Float64}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64
end

struct priorCovBiLinear
    A::Matrix{Float64}
    B::Array{Float64,3}
    C::Matrix{Float64}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64
end

struct priorCovNonLinear
    A::Matrix{Float64}
    B::Array{Float64,3}
    C::Matrix{Float64}
    D::Array{Float64,3}
    transit::Vector{Float64}
    decay::Vector{Float64}
    epsilon::Float64
end

#Base.@kwdef mutable struct PriorDCM
struct PriorDCMLinear
    pE::priorMeanLinear
    pC::priorCovLinear
end

struct PriorDCMBiLinear
    pE::priorMeanBiLinear
    pC::priorCovBiLinear
end

struct PriorDCMNonLinear
    pE::priorMeanNonLinear
    pC::priorCovNonLinear
end

"""
$(TYPEDEF)

Representation of an rDCM with fixed network architecture.

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct RigidRdcm <: RDCM
    "Binary indicator matrix for endogenous connectivity"
    a::BitMatrix
    "Binary indicator matrix for driving inputs"
    c::BitMatrix
    "Number of data points per region"
    scans::Int
    "Number of regions in network"
    const nr::Int
    "Input structure with information about driving input"
    U::InputU
    "Data structure containing BOLD signal"
    Y::BoldY
    "Connectivity parameters (A and C matrix)"
    Ep::TrueParamLinear #TODO: don't need this for model inversion
    "Confound structure"
    Conf::Confound
    "Hemodynamic response function"
    hrf::Vector{Float64}

    @doc """
    $(SIGNATURES)

    Constructor with sanity checks for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
    """
    function RigidRdcm(a, c, scans, nr, U, Y, Ep, Conf, hrf)
        if size(a, 1) ≠ size(c, 1) ||
            nr ≠ size(a, 1) ||
            (!isnothing(Y.y) && size(Matrix{Float64}(Y.y), 2) ≠ nr) ||
            (!isnothing(Y.y) && size(Matrix{Float64}(Y.y), 1) ≠ scans) ||
            size(Ep.A, 1) ≠ nr
            error("Dimension mismatch.")
        end

        if size(U.u, 2) ≠ size(c, 2)
            error("Number of inputs don't match.")
        end
        y = Y.y # for JET not to throw error
        if !isnothing(y)
            r_dt = 1
            try
                r_dt = Int(Y.dt / U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
            if size(y, 1) ≠ size(U.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsistent.")
            end
        end

        return new(a, c, scans, nr, U, Y, Ep, Conf, hrf)
    end
end

"""
$(TYPEDEF)

Representation of a sparse rDCM (sparsity constraints on network architecture).

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SparseRdcm <: RDCM
    "Binary indicator matrix for endogenous connectivity"
    a::BitMatrix
    "Binary indicator matrix for driving inputs"
    c::BitMatrix
    "Number of data points per region"
    scans::Int
    "Number of regions in network"
    const nr::Int
    "Input structure with information about driving input"
    U::InputU
    "Data structure containing BOLD signal"
    Y::BoldY
    "Connectivity parameters (A and C matrix)"
    Ep::TrueParamLinear
    "Confound structure"
    Conf::Confound
    "Hemodynamic response function"
    hrf::Vector{Float64}
    "Inform region specific sparseness (e.g., by anatomical information)"
    inform_p0::Bool
    "Prior belief about network sparseness"
    p0::Float64

    @doc """
    $(SIGNATURES)

    Constructor with sanity checks for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
    """
    function SparseRdcm(a, c, scans, nr, U, Y, Ep, Conf, hrf, inform_p0, p0)
        if size(a, 1) ≠ size(c, 1) ||
            nr ≠ size(a, 1) ||
            (!isnothing(Y.y) && size(Matrix{Float64}(Y.y), 2) ≠ nr) ||
            (!isnothing(Y.y) && size(Matrix{Float64}(Y.y), 1) ≠ scans) ||
            size(Ep.A, 1) ≠ nr
            error("Dimension mismatch.")
        end

        if size(U.u, 2) ≠ size(c, 2)
            error("Number of inputs don't match.")
        end

        if p0 < 0.0 || p0 > 1.0
            error("p0 is not a proper Bernoulli parameter.")
        end

        y = Y.y # for JET not to throw error
        if !isnothing(y)
            r_dt = 1
            try
                r_dt = Int(Y.dt / U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
            if size(y, 1) ≠ size(U.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsistent.")
            end
        end
        return new(a, c, scans, nr, U, Y, Ep, Conf, hrf, inform_p0, p0)
    end
end
