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
    maxIter::Int64
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
    maxIter::Int64
    "Tolerance for convergence"
    tol::Float64
    "Number of reruns"
    reruns::Int64
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
    #filterStr::Int64 = 0
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
    verbose::Int64
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

Ouput after inversion of [`RigidRdcm`](@ref RigidRdcm).

# Fields

$(TYPEDFIELDS)
"""
struct RigidOutput <: ModelOutput
    "Negative free energy of the whole model"
    F::Float64
    "Region specific negative free energy"
    F_r::Vector{Float64}
    "Number of iterations per region until convergence"
    iter_all::Vector{Int64}
    "Posterior shape parameter for noise"
    a_all::Vector{Float64}
    "Posterior rate parameter for noise"
    b_all::Vector{Float64}
    "Posterior mean for connectivity parameters"
    m_all::Matrix{Float64}
    "Posterior covariance for connectivity parameters"
    Σ_all::Vector{SparseMatrixCSC{Float64,Int64}}
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
            error("Inconsisten number of regions.")
        end
        #if any((!).(isposdef.(Σ_all)))
        #    error("One or more covariance matrices are not positive semi-definite.")
        #end
        return new(F, F_r, iter_all, a_all, b_all, m_all, Σ_all, inversion)
    end
end

"""
$(TYPEDEF)

Ouput after inversion of [`SparseRdcm`](@ref SparseRdcm).

# Fields

$(TYPEDFIELDS)
"""
struct SparseOutput <: ModelOutput
    "Negative free energy of the whole model"
    F::Float64
    "Region specific negative free energy"
    F_r::Vector{Float64}
    "Number of iterations per region until convergence"
    iter_all::Vector{Int64} # number of iterations per region until convergence
    "Posterior shape parameter for noise"
    a_all::Vector{Float64}
    "Posterior rate parameter for noise"
    b_all::Vector{Float64}
    "Posterior mean for connectivity parameters"
    m_all::Matrix{Float64}
    "Posterior covariance for connectivity parameters"
    Σ_all::Vector{SparseMatrixCSC{Float64,Int64}}
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
            error("Inconsisten number of regions.")
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
    scans::Int64
    "number of regions in network"
    const nr::Int64
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
                    error("Length of BOLD signal and driving input u is inconsisten.")
                end
            end
        end

        return new(a, c, scans, nr, U, Y, Ep, Conf)
    end
end

"""
$(TYPEDEF)

Represenetation of a bi-linear DCM.

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
    scans::Int64
    "number of regions in network"
    const nr::Int64
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
                error("Length of BOLD signal and driving input u is inconsisten.")
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
    scans::Int64
    "number of regions in network"
    const nr::Int64
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
                error("Length of BOLD signal and driving input u is inconsisten.")
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
    scans::Int64
    "Number of regions in network"
    const nr::Int64
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
                error("Length of BOLD signal and driving input u is inconsisten.")
            end
        end

        return new(a, c, scans, nr, U, Y, Ep, Conf, hrf)
    end
end

Base.@kwdef mutable struct BiLinearRigidRdcm <: RDCM
    "Binary indicator matrix for endogenous connectivity"
    a::BitMatrix
    "Binary indicator matrix for bi-linear dynamics."
    b::BitArray{3}
    "Binary indicator matrix for driving inputs"
    c::BitMatrix
    "Number of data points per region"
    scans::Int64
    "Number of regions in network"
    const nr::Int64
    "Input structure with information about driving input"
    U::InputU
    "Data structure containing BOLD signal"
    Y::BoldY
    "Connectivity parameters (A and C matrix)"
    Ep::TrueParamBiLinear #TODO: don't need this for model inversion
    "Confound structure"
    Conf::Confound
    "Hemodynamic response function"
    hrf::Vector{Float64}

    @doc """
    $(SIGNATURES)

    Constructor with sanity checks for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
    """
    function BiLinearRigidRdcm(a, b, c, scans, nr, U, Y, Ep, Conf, hrf)
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
                error("Length of BOLD signal and driving input u is inconsisten.")
            end
        end

        return new(a, b, c, scans, nr, U, Y, Ep, Conf, hrf)
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
    scans::Int64
    "Number of regions in network"
    const nr::Int64
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
                error("Length of BOLD signal and driving input u is inconsisten.")
            end
        end
        return new(a, c, scans, nr, U, Y, Ep, Conf, hrf, inform_p0, p0)
    end
end

#---------------------------------------------------------------------------
# outer constructor for RigidInversionParams
#---------------------------------------------------------------------------
"""
    RigidInversionParams(; maxIter=500, tol=1.0e-5)

Constructor for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
"""
function RigidInversionParams(; maxIter=500, tol=1.0e-5)
    return RigidInversionParams(maxIter, tol)
end

#---------------------------------------------------------------------------
# outer constructor for SparseInversionParams
#---------------------------------------------------------------------------
"""
    SparseInversionParams(; maxIter=500, tol=1.0e-5, reruns=100, restrictInputs=true)

Constructor for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
"""
function SparseInversionParams(; maxIter=500, tol=1.0e-5, reruns=100, restrictInputs=true)
    return SparseInversionParams(maxIter, tol, reruns, restrictInputs)
end

#---------------------------------------------------------------------------
# outer constructors for InputU
#---------------------------------------------------------------------------
function InputU(u::Matrix{Float64}, dt::Float64)
    k = size(u, 2)
    names = ["u_" * string(i) for i in 1:k]
    return InputU(u, dt, names)
end

function InputU(
    u::SparseArrays.SparseMatrixCSC{Float64,Int64}, dt::Float64, name::Vector{Any}
)
    return InputU(Matrix(u), dt, string.(name))
end

function InputU(u::Matrix{Float64}, dt::Float64, name::Vector{Any})
    return InputU(Matrix(u), dt, string.(name))
end

#---------------------------------------------------------------------------
# outer constructors for Confound
#---------------------------------------------------------------------------
function Confound(X0::Union{Matrix{Float64},Vector{Float64}}, name::Vector{Any})
    return Confound(X0, string.(name))
end

function Confound(X0::Union{Matrix{Float64},Vector{Float64}})
    nc = size(X0, 2)
    names = ["Conf_" * string(i) for i in 1:nc]
    return Confound(X0, names)
end

#---------------------------------------------------------------------------
# outer constructors and setter function for BoldY
#---------------------------------------------------------------------------
"""
$(TYPEDSIGNATURES)

Create data structure containing BOLD sigal `y`, sampling interval `dt` in seconds and optionaly `name` specifying the names of the regions.

# Examples
```julia-repl
julia> BoldY(zeros(100,3),0.5;name=["region1","region2","region3"])
```
"""
function BoldY(
    y::Union{Matrix{Float64},Nothing},
    dt::Float64;
    name::Union{Vector{Any},Vector{String},Nothing}=nothing,
)
    if isnothing(name)
        if !isnothing(y)
            nr = size(y, 2)
            name = ["y_" * string(i) for i in 1:nr]
        else
            return BoldY(y, dt, name)
        end
    end

    return BoldY(y, dt, string.(name))
end

function Base.setproperty!(val::BoldY, key::Symbol, x)
    if key == :y
        if isnothing(x) && !isnothing(val.name)
            setfield!(val, :name, nothing)
        end
        #if !isnothing(x) && isnothing(val.name)
        #    nr = size(x,2)
        #    setfield!(val,:name,["y_"*string(i) for i=1:nr])
        #end
        if !isnothing(x) && !isnothing(val.name)
            if size(x, 2) ≠ length(Vector{String}(val.name))
                error("Size of y and name vector don't match.")
            end
        end
        setfield!(val, :y, x)
    elseif key == :name
        if isnothing(x) && !isnothing(val.y)
            error("Cannot set name to nothing because y is not nothig.")
        end
        if !isnothing(x) && isnothing(val.y)
            error("Cannot set name property because y is nothing.")
        end
        if !isnothing(x) && !isnothing(val.y)
            if size(Matrix{Float64}(val.y), 2) ≠ length(x)
                error("Size of y and name vector don't match.")
            end
        end
        setfield!(val, :name, x)
    elseif key == :dt
        if x < 0.0
            error("Sampling rate needs to be positive.")
        end
        setfield!(val, :dt, x)
    end
end

#---------------------------------------------------------------------------
# outer constructors for TrueParamLinear
#---------------------------------------------------------------------------
function TrueParamLinear(A::Matrix{Float64}, C::Matrix{Float64})
    prior = get_prior_stats(A .≠ 0, C .≠ 0)
    transit = prior.pE.transit
    decay = prior.pE.decay
    epsilon = prior.pE.epsilon
    return TrueParamLinear(A, C, transit, decay, epsilon)
end

function TrueParamLinear(
    a::BitMatrix, c::BitMatrix; sample=true, fixHRF=true, rng=Xoshiro()
)
    prior = get_prior_stats(a, c)

    if sample
        A, C, transit, decay, epsilon = sample_from_prior(a, c, prior; fixHRF, rng=rng)
        TrueParamLinear(A, C, transit, decay, epsilon)
    else
        pE = prior.pE
        TrueParamLinear(pE.A, pE.C, pE.transit, pE.decay, pE.epsilon)
    end
end

function TrueParamLinear(tp::TrueParamBiLinear)
    return TrueParamLinear(tp.A, tp.C, tp.transit, tp.decay, tp.epsilon)
end

function TrueParamLinear(tp::TrueParamNonLinear)
    return TrueParamLinear(tp.A, tp.C, tp.transit, tp.decay, tp.epsilon)
end

#---------------------------------------------------------------------------
# outer constructors for TrueParamBiLinear
#---------------------------------------------------------------------------
function TrueParamBiLinear(A::Matrix{Float64}, B::Array{Float64}, C::Matrix{Float64})
    prior = get_prior_stats(A .≠ 0, B .≠ 0, C .≠ 0)
    transit = prior.pE.transit
    decay = prior.pE.decay
    epsilon = prior.pE.epsilon
    return TrueParamBiLinear(A, B, C, transit, decay, epsilon)
end

function TrueParamBiLinear(
    a::BitMatrix, b::BitArray{3}, c::BitMatrix; sample=true, fixHRF=true, rng=Xoshiro()
)
    prior = get_prior_stats(a, b, c)
    if sample
        A, B, C, transit, decay, epsilon = sample_from_prior(
            a, b, c, prior; fixHRF, rng=rng
        )
        TrueParamBiLinear(A, B, C, transit, decay, epsilon)
    else
        pE = prior.pE
        TrueParamBiLinear(pE.A, pE.B, pE.C, pE.transit, pE.decay, pE.epsilon)
    end
end

function TrueParamBiLinear(Ep::TrueParamLinear, B::Array{Float64,3})
    return TrueParamBiLinear(Ep.A, B, Ep.C, Ep.transit, Ep.decay, Ep.epsilon)
end

#---------------------------------------------------------------------------
# outer constructors for TrueParamNonLinear
#---------------------------------------------------------------------------
function TrueParamNonLinear(A::Matrix{Float64}, C::Matrix{Float64})
    nr = size(A, 1)
    nu = size(C, 2)
    B = zeros(nr, nr, nu)
    D = zeros(nr, nr, nr)
    prior = get_prior_stats(A .≠ 0, C .≠ 0)
    transit = prior.pE.transit
    decay = prior.pE.decay
    epsilon = prior.pE.epsilon
    return TrueParamNonLinear(A, B, C, D, transit, decay, epsilon)
end

function TrueParamNonLinear(
    a::BitMatrix,
    b::BitArray{3},
    c::BitMatrix,
    d::BitArray{3};
    sample=true,
    fixHRF=true,
    rng=Xoshiro(),
)
    prior = get_prior_stats(a, b, c, d)

    if sample
        A, B, C, D, transit, decay, epsilon = sample_from_prior(
            a, b, c, d, prior; fixHRF, rng=rng
        )
        TrueParamNonLinear(A, B, C, D, transit, decay, epsilon)
    else
        pE = prior.pE
        TrueParamNonLinear(pE.A, pE.B, pE.C, pE.D, pE.transit, pE.decay, pE.epsilon)
    end
end

function TrueParamNonLinear(Ep::TrueParamLinear)
    nr = size(Ep.A, 1)
    nu = size(Ep.C, 2)
    B = zeros(nr, nr, nu)
    D = zeros(nr, nr, nr)
    return TrueParamNonLinear(Ep.A, B, Ep.C, D, Ep.transit, Ep.decay, Ep.epsilon)
end

function TrueParamNonLinear(Ep::TrueParamBiLinear)
    nr = size(Ep.A, 1)
    D = zeros(nr, nr, nr)
    return TrueParamNonLinear(Ep.A, Ep.B, Ep.C, D, Ep.transit, Ep.decay, Ep.epsilon)
end

#---------------------------------------------------------------------------
# outer constructors and setter function for LinearDCM
#---------------------------------------------------------------------------
"""
$(TYPEDSIGNATURES)

Create linear DCM.

# Arguments
- `a::Matrix`: Binary indicator matrix for endogenous connectivity
- `c::Matrix`: Binary indicator matrix for driving inputs
- `scans::Int64`: Number of data points per region
- `nr::Int64`: Number of regions
- `U::InputU`: Input structure
- `Y::Union{BoldY,Nothing}`: Data structure containing BOLD signal, can be nothing
- `EP::TrueParamLinear`: Connectivity parameters containing A and C matrix
"""
function LinearDCM(
    a::Matrix{T1},
    c::Matrix{T2},
    scans::Int64,
    nr::Int64,
    U::Union{InputU,Nothing},
    Y::Union{BoldY,Nothing},
    Ep::TrueParamLinear,
) where {T1<:Number,T2<:Number}
    return LinearDCM(a .≠ 0, c .≠ 0, scans, nr, U, Y, Ep, nothing)
end

function LinearDCM(
    a::T1,
    c::T2,
    scans::Int64,
    nr::Int64,
    U::Union{InputU,Nothing},
    Y::Union{BoldY,Nothing},
    Ep::TrueParamLinear,
) where {T1<:Number,T2<:Number}
    return LinearDCM(
        BitArray(zeros(1, 1) .+ abs(a)),
        BitArray(zeros(1, 1) .+ abs(c)),
        scans,
        nr,
        U,
        Y,
        Ep,
        nothing,
    )
end

function LinearDCM(
    a::BitMatrix,
    c::BitMatrix,
    scans::Int64,
    nr::Int64,
    U::Union{InputU,Nothing},
    Y::Union{BoldY,Nothing},
    Ep::TrueParamLinear,
)
    return LinearDCM(a, c, scans, nr, U, Y, Ep, nothing)
end

"""
$(TYPEDSIGNATURES)

Create linear DCM from a bilinear or nonlinear DCM `dcm`.
"""
function LinearDCM(dcm::T) where {T<:DCM}
    return LinearDCM(
        dcm.a, dcm.c, dcm.scans, dcm.nr, dcm.U, dcm.Y, TrueParamLinear(dcm.Ep), dcm.Conf
    )
end

"""
$(TYPEDSIGNATURES)

Create linear DCM from an rDCM structure `rdcm` and output from model inversion `output`.
"""
function LinearDCM(rdcm::T1, output::T2) where {T1<:RDCM} where {T2<:ModelOutput}
    if rdcm isa RigidRdcm && output isa SparseOutput
        error("The output is from a sparse rDCM but the rDCM struct is a rigid rDCM.")
    end
    if rdcm isa SparseRdcm && output isa RigidOutput
        error("The output is from a rigid rDCM but the rDCM struct is a sparse rDCM.")
    end
    Y = BoldY(nothing, rdcm.Y.dt, rdcm.Y.name)
    nr = rdcm.nr
    A = output.m_all[:, 1:nr]
    C = output.m_all[:, (nr + 1):end]
    Ep = TrueParamLinear(A, C, rdcm.Ep.transit, rdcm.Ep.decay, rdcm.Ep.epsilon)
    return LinearDCM(rdcm.a, rdcm.c, rdcm.scans, rdcm.nr, rdcm.U, Y, Ep, rdcm.Conf)
end

function Base.setproperty!(val::LinearDCM, key::Symbol, x)
    if key == :a
        if size(x, 1) ≠ val.nr || size(x, 1) ≠ size(x, 2)
            error("Number of regions does not match.")
        end
        setfield!(val, :a, x)
    elseif key == :c
        if size(x, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x, 2) ≠ size(val.U.u, 2)
            error("Number of inputs does not match.")
        end
        setfield!(val, :c, x)
    elseif key == :scans
        if x ≠ size(Matrix{Float64}(val.Y.y), 1)
            error("Number of scans does not match.")
        end
        setfield!(val, :scans, x)
    elseif key == :U
        if !isnothing(x)
            if size(x.u, 2) ≠ size(val.c, 2)
                error("Number of inputs does not match.")
            end
            r_dt = 1
            try
                r_dt = Int(val.Y.dt / x.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
            y = val.Y.y
            if !isnothing(y)
                if size(y, 1) ≠ size(x.u, 1) / r_dt
                    error("Length of BOLD signal and driving input u is inconsisten.")
                end
            end
        end
        setfield!(val, :U, x)
    elseif key == :Y
        if !isnothing(x) && !isnothing(x.y)
            if size(x.y, 2) ≠ val.nr
                error("Number of regions does not match.")
            end
            if size(x.y, 1) ≠ val.scans
                error("Number of scans does not match.")
            end
            r_dt = 1
            try
                r_dt = Int(x.dt / val.U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
        end
        setfield!(val, :Y, x)
    elseif key == :Ep
        if size(x.A, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x.C, 2) ≠ size(val.c, 2)
            error("Number of inputs does not match.") #TODO: create macros for this error messages
        end
        setfield!(val, :Ep, x)
    end
end

#---------------------------------------------------------------------------
# outer constructors and setter function for BiLinearDCM
#---------------------------------------------------------------------------
function BiLinearDCM(rdcm::BiLinearRigidRdcm, output::T) where {T<:ModelOutput}
    Y = BoldY(nothing, rdcm.Y.dt, rdcm.Y.name)
    nr = rdcm.nr
    nu = size(rdcm.U.u,2)
    A = output.m_all[:, 1:nr]
    B = reshape(output.m_all[:, (nr + 1):(nu+1)*nr],nr,nr,nu)
    C = output.m_all[:, (nu + 1)*nr+1:end]
    Ep = TrueParamBiLinear(A, B, C, rdcm.Ep.transit, rdcm.Ep.decay, rdcm.Ep.epsilon)
    return BiLinearDCM(rdcm.a, rdcm.b, rdcm.c, rdcm.scans, rdcm.nr, rdcm.U, Y, Ep, rdcm.Conf)
end

function BiLinearDCM(dcm::LinearDCM)
    nu = size(dcm.c, 2)
    nr = size(dcm.a, 1)
    B = zeros(nr, nr, nu)
    if !isnothing(dcm.Y)
        dcm.Y.y = nothing
    end
    return BiLinearDCM(
        dcm.a,
        B .≠ 0.0,
        dcm.c,
        dcm.scans,
        dcm.nr,
        dcm.U,
        dcm.Y,
        TrueParamBiLinear(dcm.Ep, B),
        dcm.Conf,
    )
end

function BiLinearDCM(dcm::LinearDCM, B::Array{Float64,3})
    if !isnothing(dcm.Y)
        dcm.Y.y = nothing
    end
    return BiLinearDCM(
        dcm.a,
        B .≠ 0.0,
        dcm.c,
        dcm.scans,
        dcm.nr,
        dcm.U,
        dcm.Y,
        TrueParamBiLinear(dcm.Ep, B),
        dcm.Conf,
    )
end

function Base.setproperty!(val::BiLinearDCM, key::Symbol, x)
    if key == :a
        if size(x, 1) ≠ val.nr || size(x, 1) ≠ size(x, 2)
            error("Number of regions does not match.")
        end
        setfield!(val, :a, x)
    elseif key == :b
        if size(x, 1) ≠ val.nr || size(x, 1) ≠ size(x, 2)
            error("Number of regions does not match.")
        end
        if size(x, 3) ≠ size(val.c, 2) || size(x, 3) ≠ size(val.U.u, 2)
            error("Number of inputs does not match.")
        end
        setfield!(val, :b, x)
    elseif key == :c
        if size(x, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x, 2) ≠ size(val.U.u, 2) || size(x, 2) ≠ size(val.b, 3)
            error("Number of inputs does not match.")
        end
        setfield!(val, :c, x)
    elseif key == :scans
        if x ≠ size(Matrix{Float64}(val.Y.y), 1)
            error("Number of scans does not match.")
        end
        setfield!(val, :scans, x)
    elseif key == :U
        if size(x.u, 2) ≠ size(val.c, 2)
            error("Number of inputs does not match.")
        end
        r_dt = 1
        try
            r_dt = Int(val.Y.dt / x.dt)
        catch
            error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
            of the input U (u_dt). Cannot proceed.")
        end
        y = val.Y.y
        if !isnothing(y)
            if size(y, 1) ≠ size(x.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsisten.")
            end
        end
        setfield!(val, :U, x)
    elseif key == :Y
        if !isnothing(x) && !isnothing(x.y)
            if size(x.y, 2) ≠ val.nr
                error("Number of regions does not match.")
            end
            if size(x.y, 1) ≠ val.scans
                error("Number of scans does not match.")
            end
            r_dt = 1
            try
                r_dt = Int(x.dt / val.U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
        end
        setfield!(val, :Y, x)
    elseif key == :Ep
        if size(x.A, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x.C, 2) ≠ size(val.c, 2)
            error("Number of inputs does not match.") #TODO: create macros for this error messages
        end
        setfield!(val, :Ep, x)
    end
end

#---------------------------------------------------------------------------
# outer constructors and setter function for NonLinearDCM
#---------------------------------------------------------------------------
function NonLinearDCM(dcm::LinearDCM)
    nu = size(dcm.c, 2)
    nr = size(dcm.a, 1)
    b = BitArray(zeros(nr, nr, nu))
    d = BitArray(zeros(nr, nr, nr))
    if !isnothing(dcm.Y)
        dcm.Y.y = nothing
    end
    return NonLinearDCM(
        dcm.a,
        b,
        dcm.c,
        d,
        dcm.scans,
        dcm.nr,
        dcm.U,
        dcm.Y,
        TrueParamNonLinear(dcm.Ep),
        dcm.Conf,
    )
end

function NonLinearDCM(dcm::BiLinearDCM)
    nr = size(dcm.a, 1)
    d = BitArray(zeros(nr, nr, nr))
    if !isnothing(dcm.Y)
        dcm.Y.y = nothing
    end
    return NonLinearDCM(
        dcm.a,
        dcm.b,
        dcm.c,
        d,
        dcm.scans,
        dcm.nr,
        dcm.U,
        dcm.Y,
        TrueParamNonLinear(dcm.Ep),
        dcm.Conf,
    )
end

function Base.setproperty!(val::NonLinearDCM, key::Symbol, x)
    if key == :a
        if size(x, 1) ≠ val.nr || size(x, 1) ≠ size(x, 2)
            error("Number of regions does not match.")
        end
        setfield!(val, :a, x)
    elseif key == :b
        if size(x, 1) ≠ val.nr || size(x, 1) ≠ size(x, 2)
            error("Number of regions does not match.")
        end
        if size(x, 3) ≠ size(val.c, 2) || size(x, 3) ≠ size(val.U.u, 2)
            error("Number of inputs does not match.")
        end
        setfield!(val, :b, x)
    elseif key == :c
        if size(x, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x, 2) ≠ size(val.U.u, 2) || size(x, 2) ≠ size(val.b, 3)
            error("Number of inputs does not match.")
        end
        setfield!(val, :c, x)
    elseif key == :d
        if size(x, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x, 1) ≠ size(x, 2) || size(x, 2) ≠ size(x, 3)
            error("Number of regions does not match.")
        end
        setfield!(val, :d, x)
    elseif key == :scans
        if x ≠ size(Matrix{Float64}(val.Y.y), 1)
            error("Number of scans does not match.")
        end
        setfield!(val, :scans, x)
    elseif key == :U
        if size(x.u, 2) ≠ size(val.c, 2)
            error("Number of inputs does not match.")
        end
        r_dt = 1
        try
            r_dt = Int(val.Y.dt / x.dt)
        catch
            error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
            of the input U (u_dt). Cannot proceed.")
        end
        y = val.Y.y
        if !isnothing(y)
            if size(y, 1) ≠ size(x.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsisten.")
            end
        end
        setfield!(val, :U, x)
    elseif key == :Y
        if !isnothing(x) && !isnothing(x.y)
            if size(x.y, 2) ≠ val.nr
                error("Number of regions does not match.")
            end
            if size(x.y, 1) ≠ val.scans
                error("Number of scans does not match.")
            end
            r_dt = 1
            try
                r_dt = Int(x.dt / val.U.dt)
            catch
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                of the input U (u_dt). Cannot proceed.")
            end
        end
        setfield!(val, :Y, x)
    elseif key == :Ep
        if size(x.A, 1) ≠ val.nr
            error("Number of regions does not match.")
        end
        if size(x.C, 2) ≠ size(val.c, 2)
            error("Number of inputs does not match.") #TODO: create macros for this error messages
        end
        setfield!(val, :Ep, x)
    end
end

#---------------------------------------------------------------------------
# outer constructors for RigidRdcm
#---------------------------------------------------------------------------
"""
$(TYPEDSIGNATURES)

Construct a [`$(FUNCTIONNAME)`](@ref) based on a linear DCM.
"""
function RigidRdcm(dcm::LinearDCM)
    hrf = get_hrf(dcm)
    # Conf is nothing, U exists
    if isnothing(dcm.Conf) & !isnothing(dcm.U)
        Conf = Confound(ones(size(dcm.U.u, 1)), ["Constant"])
        return RigidRdcm(dcm.a, dcm.c, dcm.scans, dcm.nr, dcm.U, dcm.Y, dcm.Ep, Conf, hrf)
        # Conf exists, U exists
    elseif !isnothing(dcm.Conf) & !isnothing(dcm.U)
        return RigidRdcm(
            dcm.a, dcm.c, dcm.scans, dcm.nr, dcm.U, dcm.Y, dcm.Ep, dcm.Conf, hrf
        )
    end

    # U is nothing
    r_dt = 16 # assumes microtime resolution is 16
    y = dcm.Y.y
    if !isnothing(y) # need to do this otherwise JET.jl gives a false positive (1/2 union split)
        N = size(y, 1) * r_dt
    end
    u_dt = dcm.Y.dt / r_dt
    U = InputU(zeros(N, 0), u_dt)
    c = BitMatrix(zeros(dcm.nr, size(U.u, 2)))
    Conf = dcm.Conf
    # Conf is nothing
    if isnothing(Conf)
        Conf = Confound(zeros(N, 0))
    end
    # Conf exist
    return RigidRdcm(dcm.a, c, dcm.scans, dcm.nr, U, dcm.Y, dcm.Ep, Conf, hrf)
end

#---------------------------------------------------------------------------
# outer constructors for BiLinearRigidRdcm
#---------------------------------------------------------------------------
"""
$(TYPEDSIGNATURES)

Construct a [`$(FUNCTIONNAME)`](@ref) based on a linear DCM.
"""
function BiLinearRigidRdcm(dcm::BiLinearDCM)
    hrf = get_hrf(dcm)
    # Conf is nothing, U exists
    if isnothing(dcm.Conf) & !isnothing(dcm.U)
        Conf = Confound(ones(size(dcm.U.u, 1)), ["Constant"])
        return BiLinearRigidRdcm(dcm.a, dcm.b, dcm.c, dcm.scans, dcm.nr, dcm.U, dcm.Y, dcm.Ep, Conf, hrf)
        # Conf exists, U exists
    elseif !isnothing(dcm.Conf) & !isnothing(dcm.U)
        return BiLinearRigidRdcm(
            dcm.a, dcm.b, dcm.c, dcm.scans, dcm.nr, dcm.U, dcm.Y, dcm.Ep, dcm.Conf, hrf
        )
    end

    # U is nothing
    r_dt = 16 # assumes microtime resolution is 16
    y = dcm.Y.y
    if !isnothing(y) # need to do this otherwise JET.jl gives a false positive (1/2 union split)
        N = size(y, 1) * r_dt
    end
    u_dt = dcm.Y.dt / r_dt
    U = InputU(zeros(N, 0), u_dt)
    c = BitMatrix(zeros(dcm.nr, size(U.u, 2)))
    Conf = dcm.Conf
    # Conf is nothing
    if isnothing(Conf)
        Conf = Confound(zeros(N, 0))
    end
    # Conf exist
    return BiLinearRigidRdcm(dcm.a, dcm.b, c, dcm.scans, dcm.nr, U, dcm.Y, dcm.Ep, Conf, hrf)
end

#---------------------------------------------------------------------------
# outer constructors for SparseRdcm
#---------------------------------------------------------------------------
"""
$(TYPEDSIGNATURES)

Construct a [`$(FUNCTIONNAME)`](@ref) based on a linear DCM.

# Arguments
- `dcm`: DCM model
- `inform_p0::Bool`: Inform region specific sparseness (e.g., by anatomical information)
- `p0::Float64`: Prior belief about network sparseness
"""
function SparseRdcm(dcm::LinearDCM; inform_p0=false, p0=0.5)
    hrf = get_hrf(dcm)
    Conf = dcm.Conf
    if isnothing(Conf)
        Conf = Confound(ones(size(dcm.U.u, 1)), ["Constant"])
        #c = [dcm.c BitMatrix(ones(dcm.nr, 1))]
    end
    return SparseRdcm(
        dcm.a, dcm.c, dcm.scans, dcm.nr, dcm.U, dcm.Y, dcm.Ep, Conf, hrf, inform_p0, p0
    )
end

#---------------------------------------------------------------------------
# outer constructors for Options
#---------------------------------------------------------------------------
"""
    Options(invParams::T;synthetic::Bool,verbose=1,testing=false,rng=Xoshiro())

Constructor for [`$(FUNCTIONNAME)`](@ref). See type description for information about arguments.
"""
function Options(
    invParams::T; synthetic::Bool, verbose=1, testing=false, rng=Xoshiro()
) where {T<:AbstractInvParam}
    return Options(invParams, verbose, synthetic, testing, rng)
end

function Base.copy(dcm::LinearDCM)
    return LinearDCM([deepcopy(getfield(dcm, k)) for k in fieldnames(LinearDCM)]...)
end
function Base.copy(dcm::BiLinearDCM)
    return BiLinearDCM([deepcopy(getfield(dcm, k)) for k in fieldnames(BiLinearDCM)]...)
end
function Base.copy(dcm::NonLinearDCM)
    return NonLinearDCM([deepcopy(getfield(dcm, k)) for k in fieldnames(NonLinearDCM)]...)
end
function Base.copy(dcm::RigidRdcm)
    return RigidRdcm([deepcopy(getfield(dcm, k)) for k in fieldnames(RigidRdcm)]...)
end
function Base.copy(dcm::SparseRdcm)
    return SparseRdcm([deepcopy(getfield(dcm, k)) for k in fieldnames(SparseRdcm)]...)
end
function Base.copy(dcm::BiLinearRigidRdcm)
    return BiLinearRigidRdcm([deepcopy(getfield(dcm, k)) for k in fieldnames(BiLinearRigidRdcm)]...)
end
