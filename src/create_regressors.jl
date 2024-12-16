
"""
    reduce_zeros!(X,Y;rng)

If there are more zero-valued frequencies than informative ones, subsamples those frequencies to balance dataset.

# Arguments
- `X::Matrix`: design matrix (predictors)
- `Y::Matrix`: data
"""
function reduce_zeros!(
    X::Matrix{T1}, Y::Matrix{T2}, rng::AbstractRNG
) where {T1<:Number,T2<:Number}
    # get all indices
    idx = collect(1:size(Y, 1))

    # data
    data = sum(abs.([Y X]); dims=2)

    # get zero frequencies
    idx_0 = idx[vec(data .== 0)]

    # number of zero frequencies
    n0 = sum(data .== 0)

    # number of non-zero frequencies
    n1 = sum(data .> 0)

    # balance the data if there are too many zeros
    if n0 > n1
        idx_del = BitVector(vec([zeros(Bool, 1, n1) ones(Bool, 1, n0 - n1)]))
        idx_del = idx_del[randperm(rng, n0)]
        Y[idx_0[idx_del], :] .= NaN
    end
    return nothing
end

function create_regressors_core(
    hrf::Vector{Float64},
    y::Matrix{Float64},
    u::Matrix{Float64},
    u_dt::Float64,
    y_dt::Float64,
    X0::Union{Vector{Float64},Matrix{Float64}},
    rng::AbstractRNG,
)
    r_dt = 1
    try
        r_dt = Int64(y_dt / u_dt)
    catch
        error(
            "The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.",
        )
    end

    Nu, nu = size(u) #input length and number of inputs
    Ny, nr = size(y)

    # Fourier transform of hemodynamic response function (HRF)
    h_fft = rfft(hrf, 1)
    # convolution of stimulus input and HRF
    uh = irfft(rfft(u, 1) .* repeat(h_fft, 1, nu), Nu, 1)

    # subsample input and HRF
    uh = uh[1:r_dt:end, :]
    X0 = X0[1:r_dt:end, :]

    YU = zeros(Ny,nu*nr)
    for j in 1:nu
        for i in 1:Ny
            YU[i,(j-1)*nr+1:j*nr] = uh[i,j] .* y[i,:]
        end
    end

    uh = [uh X0]

    # combine regressors
    X = [y[1:(end - 1), :] (YU[1:(end - 1), :] ./ r_dt) (uh[1:(end - 1), :] ./ r_dt)]
    Y = (y[2:end, :] .- y[1:(end - 1), :]) ./ y_dt

    return X, Y
end

function create_regressors!(dcm::T, rng::AbstractRNG) where {T<:RDCM}

    # if isnothing(dcm.Conf)
    #     # add only constant confound
    #     dcm.Conf = Confound(ones(size(dcm.U.u,1)),["Constant"]) # TODO: add case for resting state fMRI
    # else
    #     # add constant confound to existing confounds
    #     dcm.Conf.X0 = [dcm.Conf.X0 ones(size(dcm.U.u,1))]
    #     push!(dcm.Conf.name,"Constant")
    #     #dcm.Conf.name = [dcm.Conf.name, ["Constant"]] #TODO: add warning if there is already a constant regressor in Conf
    # end

    # # add confound regressor dimensions
    # nc = size(dcm.Conf.X0,2)
    # dcm.c = [dcm.c BitMatrix(ones(dcm.nr,nc))]

    return create_regressors_core(
        dcm.hrf, Matrix{Float64}(dcm.Y.y), dcm.U.u, dcm.U.dt, dcm.Y.dt, dcm.Conf.X0, rng
    )
end
