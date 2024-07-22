
#------------------------------------------------------------------------------------------------------------
# version which only looks at positive frequencies and also splits data into real and imaginary part
# this is the most efficient one and also correct one
# -----------------------------------------------------------------------------------------------------------

function rigid_inversion(
    rdcm::RigidRdcm, X_c::Matrix{ComplexF64}, Y_c::Matrix{ComplexF64}, opt::Options
)
    dcm = copy(rdcm) #TODO: maybe there is more elegant way where no copy is needed
    maxIter = opt.invParams.maxIter
    pr = opt.invParams.tol^2

    # TODO: might make sense to put this part into create regressors
    # split data into real and imaginary part
    if iseven(size(Y_c, 1))
        # don't need the imaginary part of the constant and nyquist frequency because it's zero
        Y = [real(Y_c); imag(Y_c)[2:(end - 1), :]]
        X = [real(X_c); imag(X_c)[2:(end - 1), :]]
    else
        # if size of Y is odd there is no nyquist frequency, so only cut away the imaginary part of const. frequency
        Y = [real(Y_c); imag(Y_c)[2:end, :]]
        X = [real(X_c); imag(X_c)[2:end, :]]
    end

    nr = size(Y, 2)
    #Ny, nr = size(Y)

    #const_freq = X_c[1,1:nr]

    # no baseline regressor for simulations, TODO: put this also in create regressor function
    if opt.synthetic
        nc = size(dcm.Conf.X0, 2)
        dcm.c[:, (end - nc + 1):end] .= false
    end

    idx = [dcm.a dcm.c]

    μ0, l0, a0, β0 = get_priors(dcm)

    # allocate memory
    F_all = zeros(Float64, nr)
    a_all = zeros(Float64, nr)
    b_all = zeros(Float64, nr)
    m_all = zeros(Float64, size(idx)) #one could maybe also make this sparse
    Σ_all = [spzeros(Float64, (size(idx, 2), size(idx, 2))) for _ in 1:nr]
    iter_all = zeros(Int64, nr)

    # array for storing predicted derivative of signal (in frequency domain)
    # yd_fft_pred = zeros(Float64, size(X,1),nr)

    prog = Progress(nr; enabled=!opt.testing)
    for r in 1:nr
        idx_y = .!isnan.(Y[:, r])
        idx_r = idx[r, :]

        N_eff = sum(idx_y)

        # remove unnecessary dimensions
        X_r = X[idx_y, idx_r]
        Y_r = Y[idx_y, r]

        # effective dimensionality
        dim_r = sum(idx_r)

        # prior precision matrix
        l0_r = diagm(l0[r, idx_r])

        # prior mean for connectivity
        μ0_r = μ0[r, idx_r]

        # precompute X'X and X'Y
        W = X_r' * X_r
        V = X_r' * Y_r

        # initialise posterior mean of gamma distribution
        τ = a0 / β0

        # calculate posterior shape parameter
        a_r = a0 + N_eff * 0.5

        # set old F
        F_old = -Inf

        # define variable
        β_r = 0.0
        μ_r = zeros(size(μ0_r))
        F_r = 0.0
        Σ_r = zeros(size(l0_r))

        for i in 1:maxIter
            β_r, QF, τ = update_posterior_rigid!(
                μ_r, Σ_r, a_r, τ, W, l0_r, μ0_r, V, Y_r, X_r, β0
            )

            F_r = compute_F(N_eff, a_r, β_r, QF, τ, l0_r, μ_r, μ0_r, Σ_r, a0, β0, dim_r)

            # check for convergence
            if (F_old - F_r)^2 < pr
                iter_all[r] = i
                break
            end

            # store old negatve free energy
            F_old = F_r
        end

        # store region specific negative free energy
        F_all[r] = F_r
        # store region specific parameter estimates
        a_all[r] = a_r
        b_all[r] = β_r
        m_all[r, idx_r] = μ_r
        Σ_all[r][idx_r, idx_r] = Σ_r

        # get the predicted signal from the GLM (in frequency domain)
        # yd_fft_pred[:,r] = X[:,idx_r] * μ_r

        # update progressbar
        next!(prog)
    end

    nc = size(dcm.Conf.X0, 2)
    m_all = m_all[:, 1:(end - nc)] # cut away regressor estimate
    for r in 1:nr
        Σ_all[r] = Σ_all[r][1:(end - nc), 1:(end - nc)] # TODO: make a test for this
    end

    # y_pred = convert_to_time_domain(yd_fft_pred,Ny,dcm.Y.dt,const_freq) it's simply not possible to get good predictions from predicting the derivative in frequency domain

    return RigidOutput(
        sum(F_all), F_all, iter_all, a_all, b_all, m_all, Σ_all, "tapas_rdcm_rigid"
    )
end

function update_posterior_rigid!(
    μ_r::Vector{Float64},
    Σ_r::Matrix{Float64},
    a_r::Float64,
    τ::Float64,
    W::Matrix{Float64},
    l0_r::Matrix{Float64},
    μ0_r::Vector{Float64},
    V::Vector{Float64},
    Y_r::Vector{Float64},
    X_r::Matrix{Float64},
    β0::Float64,
)

    # update posterior covariance matrix
    Σ_r .= inv(τ * W + l0_r)

    # update posterior mean
    μ_r .= Σ_r * (τ * V + l0_r * μ0_r)

    # update posterior rate parameter
    QF = ((Y_r - X_r * μ_r)' * (Y_r - X_r * μ_r) + tr(X_r' * X_r * Σ_r)) * 0.5
    #QF = (Y_r'*Y_r - 2μ_r'*V + μ_r'*W*μ_r + tr(W*Σ_r))*0.5 # TODO: this line is probably more efficient
    β_r = β0 + QF

    # update posterior mean of Gamma distribution
    τ = a_r / β_r
    return β_r, QF, τ
end

function compute_F(
    N_eff::Int64,
    a_r::Float64,
    β_r::Float64,
    QF::Float64,
    τ::Float64,
    l0_r::Matrix{Float64},
    μ_r::Vector{Float64},
    μ0_r::Vector{Float64},
    Σ_r::Matrix{Float64},
    a0::Float64,
    β0::Float64,
    dim_r::Int64,
)

    # compute components of negative free energy
    log_lik = 0.5 * (N_eff * (digamma(a_r) - log(β_r)) - N_eff * log(2π)) - QF * τ
    log_p_weight =
        0.5 * (
            logdet(l0_r) - dim_r * log(2π) - (μ_r - μ0_r)' * l0_r * (μ_r - μ0_r) -
            tr(l0_r * Σ_r)
        ) #TODO: check if logdet is the most efficient way and also doesn't give the same results for not positive define matrices as matlab
    log_p_prec =
        a0 * log(β0) - loggamma(a0) + (a0 - 1.0) * (digamma(a_r) - log(β_r)) - β0 * τ
    log_q_weight = 0.5 * (logdet(Σ_r) + dim_r * (1.0 + log(2π)))
    log_q_prec = a_r - log(β_r) + loggamma(a_r) + (1.0 - a_r) * digamma(a_r)

    # compute region specific negative free energy
    return log_lik + log_p_weight + log_p_prec + log_q_weight + log_q_prec
end
