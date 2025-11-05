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
    u::SparseArrays.SparseMatrixCSC{Float64,Int}, dt::Float64, name::Vector{Any}
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
            error("Cannot set name to nothing because y is not nothing.")
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
- `scans::Int`: Number of data points per region
- `nr::Int`: Number of regions
- `U::InputU`: Input structure
- `Y::Union{BoldY,Nothing}`: Data structure containing BOLD signal, can be nothing
- `EP::TrueParamLinear`: Connectivity parameters containing A and C matrix
"""
function LinearDCM(
    a::Matrix{T1},
    c::Matrix{T2},
    scans::Int,
    nr::Int,
    U::Union{InputU,Nothing},
    Y::Union{BoldY,Nothing},
    Ep::TrueParamLinear,
) where {T1<:Number,T2<:Number}
    return LinearDCM(a .≠ 0, c .≠ 0, scans, nr, U, Y, Ep, nothing)
end

function LinearDCM(
    a::T1,
    c::T2,
    scans::Int,
    nr::Int,
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
    scans::Int,
    nr::Int,
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
    A = output.μ[:, 1:nr]
    C = output.μ[:, (nr + 1):end]
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
                #! format: off
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
                #! format: on
            end
            y = val.Y.y
            if !isnothing(y)
                if size(y, 1) ≠ size(x.u, 1) / r_dt
                    error("Length of BOLD signal and driving input u is inconsistent.")
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
                #! format: off
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
                #! format: on
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
            #! format: off
            error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
            #! format: on
        end
        y = val.Y.y
        if !isnothing(y)
            if size(y, 1) ≠ size(x.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsistent.")
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
                #! format: off
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
                #! format: on
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
            #! format: off
            error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
            #! format: on
        end
        y = val.Y.y
        if !isnothing(y)
            if size(y, 1) ≠ size(x.u, 1) / r_dt
                error("Length of BOLD signal and driving input u is inconsistent.")
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
                #! format: off
                error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
                #! format: on
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

#---------------------------------------------------------------------------
# Implementations of copy for custom structs
#---------------------------------------------------------------------------

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

"""
    isapproxsymmetric(A)

Check symmetry of matrix A.
"""
function isapproxsymmetric(A)
    return isapprox(A, A'; rtol=1e-12)
end
