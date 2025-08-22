
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
    nc = size(dcm.Conf.X0, 2)
    if nc == 1
        conf_weight_idx = BitVector(ones(nr))
    else
        conf_weight_idx = BitMatrix(ones(nr, nc))
    end
    #if opt.synthetic
    #    dcm.c[:, (end - nc + 1):end] .= false
    #end

    idx = [dcm.a dcm.c conf_weight_idx]

    μ0, l0, a0, β0 = get_priors(dcm, conf_weight_idx)

    # allocate memory
    F_all = zeros(Float64, nr)
    a_all = zeros(Float64, nr)
    b_all = zeros(Float64, nr)
    m_all = zeros(Float64, size(idx)) #one could maybe also make this sparse
    Σ_all = [spzeros(Float64, (size(idx, 2), size(idx, 2))) for _ in 1:nr]
    iter_all = ones(Int, nr)

    prog = Progress(nr; enabled=(!opt.testing))
    for r in 1:nr
        idx_y = .!isnan.(Y[:, r])
        idx_r = idx[r, :]

        N_eff = sum(idx_y)

        # remove unnecessary dimensions
        Xᵣ = X[idx_y, idx_r]
        Yᵣ = Y[idx_y, r]

        # effective dimensionality
        dim_r = sum(idx_r)

        # prior precision matrix
        l0ᵣ = diagm(l0[r, idx_r])

        # prior mean for connectivity
        μ0ᵣ = μ0[r, idx_r]

        # precompute X'X and X'Y
        W = Xᵣ' * Xᵣ
        V = Xᵣ' * Yᵣ

        # initialise posterior mean of gamma distribution
        τ = a0 / β0

        # calculate posterior shape parameter
        aᵣ = a0 + N_eff * 0.5

        # set old F
        F_old = -Inf

        # define variable
        βᵣ = 0.0
        μᵣ = zeros(size(μ0ᵣ))
        Fᵣ = 0.0
        Σᵣ = zeros(size(l0ᵣ))

        for i in 1:maxIter
            βᵣ, QF, τ = update_posterior_rigid!(μᵣ, Σᵣ, aᵣ, τ, W, l0ᵣ, μ0ᵣ, V, Yᵣ, Xᵣ, β0)

            Fᵣ = compute_F(N_eff, aᵣ, βᵣ, QF, τ, l0ᵣ, μᵣ, μ0ᵣ, Σᵣ, a0, β0, dim_r)

            # check for convergence
            if (F_old - Fᵣ)^2 < pr
                iter_all[r] = i
                break
            end

            if i == maxIter
                @warn "Reached maximum number of iterations for region $(r)."
            end

            # store old negative free energy
            F_old = Fᵣ
        end

        # store region specific negative free energy
        F_all[r] = Fᵣ
        # store region specific parameter estimates
        a_all[r] = aᵣ
        b_all[r] = βᵣ
        m_all[r, idx_r] = μᵣ
        Σ_all[r][idx_r, idx_r] = Σᵣ

        # update progressbar
        next!(prog)
    end

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
    μᵣ::Vector{Float64},
    Σᵣ::Matrix{Float64},
    a_r::Float64,
    τᵣ::Float64,
    W::Matrix{Float64},
    l0ᵣ::Matrix{Float64},
    μ0ᵣ::Vector{Float64},
    V::Vector{Float64},
    Yᵣ::Vector{Float64},
    Xᵣ::Matrix{Float64},
    β0::Float64,
)

    # update posterior covariance matrix
    Σᵣ .= inv(τᵣ * W + l0ᵣ)

    # update posterior mean
    μᵣ .= Σᵣ * (τᵣ * V + l0ᵣ * μ0ᵣ)

    # update posterior rate parameter
    QF = ((Yᵣ - Xᵣ * μᵣ)' * (Yᵣ - Xᵣ * μᵣ) + tr(Xᵣ' * Xᵣ * Σᵣ)) * 0.5
    #QF = (Yᵣ'*Yᵣ - 2μᵣ'*V + μᵣ'*W*μᵣ + tr(W*Σᵣ))*0.5 # TODO: this line is probably more efficient
    βᵣ = β0 + QF

    # update posterior mean of Gamma distribution
    τᵣ = a_r / βᵣ
    return βᵣ, QF, τᵣ
end

function compute_F(
    N_eff::Int,
    aᵣ::Float64,
    βᵣ::Float64,
    QF::Float64,
    τᵣ::Float64,
    l0ᵣ::Matrix{Float64},
    μᵣ::Vector{Float64},
    μ0ᵣ::Vector{Float64},
    Σᵣ::Matrix{Float64},
    a0::Float64,
    β0::Float64,
    dimᵣ::Int,
)

    # compute components of negative free energy
    log_lik = 0.5 * (N_eff * (digamma(aᵣ) - log(βᵣ)) - N_eff * log(2π)) - QF * τᵣ
    log_p_weight =
        0.5 * (logdet(l0ᵣ) - dimᵣ * log(2π) - (μᵣ - μ0ᵣ)' * l0ᵣ * (μᵣ - μ0ᵣ) - tr(l0ᵣ * Σᵣ)) #TODO: check if logdet is the most efficient way and also doesn't give the same results for not positive define matrices as matlab
    log_p_prec =
        a0 * log(β0) - loggamma(a0) + (a0 - 1.0) * (digamma(aᵣ) - log(βᵣ)) - β0 * τᵣ
    log_q_weight = 0.5 * (logdet(Σᵣ) + dimᵣ * (1.0 + log(2π)))
    log_q_prec = aᵣ - log(βᵣ) + loggamma(aᵣ) + (1.0 - aᵣ) * digamma(aᵣ)

    # compute region specific negative free energy
    return log_lik + log_p_weight + log_p_prec + log_q_weight + log_q_prec
end
