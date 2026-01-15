"""
    generate_BOLD(dcm; SNR, TR=NaN, rng=Xoshiro(), triple_input=true)

Generate synthetic BOLD signal timeseries based on a DCM.

# Arguments
- `dcm::T where T <: DCM`: DCM structure (can be a linear, bilinear on nonlinear DCM)
- `SNR::Real`: Signal to noise ratio
- `TR::Real`: Sampling interval in seconds (can be omitted if dt is specified in dcm.Y)
- `rng::AbstractRNG`: Random number generator for noise sampling.

# Output
- `y_noise::Matrix{Float64}`: BOLD signal timeseries with noise
- `y::Matrix{Float64}`: BOLD signal timeseries without noise
- `x::Matrix{Float64}`: Neuronal signal
- `h::Vector{Float64}`: Hemodynamic response function

# Examples
```julia-repl
julia> y_noise, y, x, h = generate_BOLD(load_example_DCM();SNR=10)
```
"""
function generate_BOLD(dcm::T; SNR::Real, TR::Real=NaN, rng=Xoshiro()) where {T<:DCM}
    dcm_c = copy(dcm)
    r_dt = 1
    if isnan(TR)
        if dcm_c.Y isa Nothing
            #! format: off
            error("Y field is empty. Please specify the TR of the BOLD signal manually using tapas_rdcm_generate(dcm;TR=value)")
            #! format: on
        else
            r_dt = Int64(dcm_c.Y.dt / dcm_c.U.dt)
        end
    else
        try
            r_dt = Int64(TR / dcm_c.U.dt)
            if dcm_c.Y isa Nothing
                dcm_c.Y = BoldY(nothing, TR)
            else
                dcm_c.Y.dt = TR
                @warn "Overwriting sampling time of Y with TR."
            end
        catch
            #! format: off
            error("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.")
            #! format: on
        end
    end

    N = size(dcm_c.U.u, 1)
    nr = size(dcm_c.a, 1)

    # allocate memory for data
    y = zeros(3*N, nr)

    # compute fixed hemodynamic response function
    h = get_hrf(N, dcm_c.U.dt)

    # compute neuronal signal x
    _, x = dcm_euler_gen(dcm_c)
    x = [x; zeros(2*N, nr)]
    h = [h; zeros(2*N)]

    # convolve neuronal signal with HRF
    for i in 1:nr
        y[:, i] = real(ifft(fft(x[:, i]) .* fft(h))) # need to take real part because imaginary part is not exactly zero due to numerical reasons
    end

    x = x[1:N, :]
    h = h[1:N, :]
    y = y[1:N, :]
    # sampling
    y = y[1:r_dt:end, :]

    # check signal
    if any(abs.(y) .== Inf) || any(isnan.(y))
        @warn "BOLD signal contains Inf or NaN. Check data generating parameters or input structure."
    elseif any(abs.(y) .> 20.0)
        @warn "BOLD signal contains large values."
    elseif all(y .== 0.0)
        @warn "BOLD signal contains only zeros. Check data generating parameters or input structure."
    end

    # add noise
    noise = randn(rng, Float64, size(y)) * diagm(vec(std(y; dims=1) ./ SNR))

    return y + noise, y, x, h
end

"""
$(SIGNATURES)

Generate a hemodynamic response function (HRF) given a certain lenght N and sampling rate u_dt.

# Arguments
- `N::Int`: Lenght of HRF
- `u_dt::Float64`: Repetition time

# Output
- `hrf::Vector{Float64}`: Hemodynamic response function
"""
function get_hrf(N::Int, u_dt::Float64)
    r_dt = 1 # in Matlab implementation this is effectively always one
    U = zeros(N, 1)
    U[1:r_dt] .= 1

    U = InputU(U, u_dt)
    Y = BoldY(nothing, u_dt)
    GT = TrueParamLinear(zeros(1, 1) .- 1, zeros(1, 1) .+ 16)

    dcm_hrf = LinearDCM(1, 1, N, 1, U, Y, GT)

    y, _ = dcm_euler_gen(dcm_hrf)

    return y[1:r_dt:end, 1]
end

"""
$(SIGNATURES)

Generate a hemodynamic response function (HRF) given a DCM.

# Arguments
- `dcm::DCM`: DCM structure (can be a linear, bilinear on nonlinear DCM)

# Output
- `hrf::Vector{Float64}`: Hemodynamic response function
"""
function get_hrf(dcm::T) where {T<:DCM}
    N = 0
    u_dt = 0
    y = nothing
    if !isnothing(dcm.Y)
        y = dcm.Y.y # need to do this otherwise JET.jl gives a false positive (1/2 union split)
    end
    if !isnothing(dcm.U)
        N = size(dcm.U.u, 1)
        u_dt = dcm.U.dt
    else
        if !isnothing(y)
            N = size(y, 1) * 16 # assumes microtime resolution is 16
        else
            error("BOLD signal is empty.")
        end
        u_dt = dcm.Y.dt / 16
    end

    return get_hrf(N, u_dt)
end
