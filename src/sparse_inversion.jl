
#------------------------------------------------------------------------------------------------------------
# version which only looks at positive frequencies and also splits data into real and imaginary part
# this is the most efficient one and also correct one
# -----------------------------------------------------------------------------------------------------------

function sparse_inversion(
    rdcm::SparseRdcm, X_c::Matrix{ComplexF64}, Y_c::Matrix{ComplexF64}, opt::Options
)
    dcm = copy(rdcm) # todo: check again why we do this

    maxIter = opt.invParams.maxIter
    pr = opt.invParams.tol
    reruns = opt.invParams.reruns

    nr = size(Y_c, 2)

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

    # get the priors
    μ0, l0, a0, b0 = get_priors(dcm, conf_weight_idx)

    # dimensionality of θ
    D = size(dcm.a, 2) + size(dcm.c, 2) + size(conf_weight_idx, 2)

    # allocate memory
    F_all = zeros(Float64, nr)
    iter_all = zeros(Int, nr)
    a_all = zeros(Float64, nr)
    b_all = zeros(Float64, nr)
    m_all = zeros(Float64, nr, D) #one could maybe also make this sparse
    Σ_all = [spzeros(Float64, (D, D)) for _ in 1:nr]
    z_all = zeros(Float64, nr, D)

    # iterate over regions
    prog = Progress(nr; enabled=(!opt.testing))
    for r in 1:nr

        # allocate memory for the results from all reruns
        μ_r_iter = zeros(D, reruns)
        Σ_r_iter = [spzeros(Float64, (D, D)) for _ in 1:reruns]
        a_r_iter = zeros(reruns)
        b_r_iter = zeros(reruns)
        iter_r_iter = ones(Int, reruns)
        # TODO: in matlab version they save every rerun also l0 (prior precision) -> this doesn't make any sense
        z_r_iter = zeros(D, reruns)
        F_r_iter = zeros(reruns)

        # remove unnecessary data points
        idx_y = .!isnan.(Y[:, r])
        X_r = X[idx_y, :] # TODO: ckeck if this can be optimised: i.e. using a view or making it immutable
        Y_r = Y[idx_y, r]

        # effective number of data points
        N_eff = sum(idx_y)

        # prior precision matrix
        l0_r = diagm(l0[r, :])

        # prior mean for connectivity
        μ0_r = μ0[r, :]

        # set p0 (Bernoulli prior)
        p0 = ones(D) * dcm.p0

        # inform p0 (e.g., by anatomical information)
        if dcm.inform_p0
            p0[1:nr] .*= dcm.a[r, :]
        end

        # ensure self-connectivity
        p0[r] = 1.0

        # ensure baseline regressor (will only contribute for empirical data)
        if opt.synthetic
            p0[end - (nc - 1)] = 0.0
        else
            p0[end - (nc - 1)] = 1.0
        end

        # make sure that driving inputs are only on correct connections
        if opt.invParams.restrictInputs
            p0[(nr + 1):(end - nc)] = dcm.c[r, :]
        end

        # allocate memory
        Σ_r = zeros(size(l0_r))
        μ_r = zeros(size(μ0_r))

        # precompute certain quantities:
        # estimate variables X'X and X'Y per region
        W = X_r' * X_r
        V = X_r' * Y_r

        for iter in 1:reruns # Matlab version starts this loop earlier but not necessary (quantities above are constant accros iterations)
            # initialise z_r, τ_r and a_r per region
            z_r = copy(p0)
            τ_r = a0 / b0
            a_r = a0 + N_eff / 2
            b_r = b0 # TODO: room for optimisation
            z_idx = BitVector(zeros(size(z_r)))
            QF = 0.0

            # define random matrix Z
            Z = diagm(z_r)

            # expectation (over Z) of ZX'XZ
            G = Z * W * Z
            G[diagind(G)] .= z_r .* diag(W)

            # set old F
            F_old = -Inf

            # convergence loop
            for i in 1:maxIter
                b_r, QF, τ_r = update_posterior_sparse!(
                    μ_r, Σ_r, a_r, τ_r, W, l0_r, μ0_r, V, Y_r, b0, z_r, p0, Z, G, D, opt
                )

                # check for sparsity (because of small values)
                μ_r[abs.(μ_r) .< 1.0e-5] .= 0.0
                z_r[μ_r .== 0.0] .= 0.0

                # get the "present" connections
                z_idx .= (z_r .> pr^2) .& (z_r .< 1.0) .> 0

                F_r = compute_F_sparse(
                    N_eff,
                    a_r,
                    b_r,
                    QF,
                    τ_r,
                    l0_r,
                    μ_r,
                    μ0_r,
                    Σ_r,
                    a0,
                    b0,
                    D,
                    z_r,
                    z_idx,
                    p0,
                )

                # check for convergence
                if (F_old - F_r)^2 < pr^2
                    iter_r_iter[iter] = i
                    break
                end

                if i == maxIter # TODO: write test
                    @warn "Reached maximum number of iterations for region $(r)."
                end

                # store old negative free energy
                F_old = F_r
            end
            # check for sparsity (because of small values)
            μ_r[abs.(μ_r) .< 1.0e-5] .= 0.0
            z_r[μ_r .== 0.0] .= 0.0

            # get the "present" connections
            z_idx .= (z_r .> pr^2) .& (z_r .< 1.0) .> 0

            F_r_iter[iter] = compute_F_sparse(
                N_eff, a_r, b_r, QF, τ_r, l0_r, μ_r, μ0_r, Σ_r, a0, b0, D, z_r, z_idx, p0
            )

            # asign iteration-specific values
            μ_r_iter[:, iter] = μ_r
            z_r_iter[:, iter] = z_r
            Σ_r_iter[iter][:, :] = Σ_r
            a_r_iter[iter] = a_r
            b_r_iter[iter] = b_r
        end

        # select the run with the highest negative free energy
        best = argmax(F_r_iter)

        # store region-specific parameters
        m_all[r, :] = μ_r_iter[:, best]
        z_all[r, :] = z_r_iter[:, best]
        Σ_all[r][:, :] = Σ_r_iter[best][:, :]
        a_all[r] = a_r_iter[best]
        b_all[r] = b_r_iter[best]
        F_all[r] = F_r_iter[best]
        iter_all[r] = iter_r_iter[best]

        # update progressbar
        next!(prog)
    end

    m_all = m_all[:, 1:(end - nc)] # cut away regressor estimate
    z_all = z_all[:, 1:(end - nc)]
    for r in 1:nr
        Σ_all[r] = Σ_all[r][1:(end - nc), 1:(end - nc)] # TODO: make a test for this
    end

    return SparseOutput(
        sum(F_all), F_all, iter_all, a_all, b_all, m_all, Σ_all, z_all, "tapas_rdcm_sparse"
    )
end

function update_posterior_sparse!(
    μ_r::Vector{Float64},
    Σ_r::Matrix{Float64},
    a_r::Float64,
    τ_r::Float64,
    W::Matrix{Float64},
    l0_r::Matrix{Float64},
    μ0_r::Vector{Float64},
    V::Vector{Float64},
    Y_r::Vector{Float64},
    b0::Float64,
    z_r::Vector{Float64},
    p0::Vector{Float64},
    Z::Matrix{Float64},
    G::Matrix{Float64},
    D::Int,
    opt::Options,
)

    # update posterior covariance matrix
    Σ_r .= inv(τ_r * G + l0_r)

    # update posterior mean
    μ_r .= Σ_r * (τ_r * Z * V + l0_r * μ0_r)

    # update posterior of binary indicator variables
    Wb = W .* (μ_r * μ_r')
    Ws = W .* Σ_r
    A = Wb + Ws

    # estimate g (as in paper)
    g = log.(p0 ./ (1.0 .- p0)) .+ τ_r * μ_r .* V .+ τ_r * diag(A) ./ 2

    if opt.testing
        # for reference testing we need to turn off permutation
        # because the RNGs of Matlab and Julia are different and it's
        # not possible to generate the same sequence of numbers
        order = 1:D
    else
        # sample a random order
        order = randperm(opt.rng, D)
    end

    # iterate over all binary variables
    for i in order
        z_r[i] = 1.0 # todo: why are we doing this?
        g[i] = g[i] - τ_r * z_r' * A[:, i]
        z_r[i] = 1.0 / (1.0 + exp(-g[i])) #this is different from matlab version because we take the real
    end

    # create Z matrix (binary indicators on diagonal)
    Z .= diagm(z_r)

    # re-estimate expectation (over Z) of ZX'XZ
    G .= Z * W * Z
    G[diagind(G)] .= z_r .* diag(W)

    # update posterior rate parameter
    QF = (Y_r' * Y_r - μ_r' * Z * V * 2.0 + μ_r' * G * μ_r + tr(G * Σ_r)) * 0.5
    b_r = b0 + QF

    # update posterior mean of Gamma distribution
    τ_r = a_r / b_r

    return b_r, QF, τ_r
end

function compute_F_sparse(
    N_eff::Int,
    a_r::Float64,
    b_r::Float64,
    QF::Float64,
    τ_r::Float64,
    l0_r::Matrix{Float64},
    μ_r::Vector{Float64},
    μ0_r::Vector{Float64},
    Σ_r::Matrix{Float64},
    a0::Float64,
    b0::Float64,
    dim_r::Int,
    z_r::Vector{Float64},
    z_idx::BitVector,
    p0::Vector{Float64},
)

    # compute components of negative free energy
    log_lik = 0.5 * (N_eff * (digamma(a_r) - log(b_r)) - N_eff * log(2π)) - τ_r * QF # TODO: In Matlab version QF is multiplied with 0.5 -> that's wrong because QF was already multiplied by 0.5
    log_p_weight =
        0.5 * (
            logdet(l0_r) - dim_r * log(2π) - (μ_r - μ0_r)' * l0_r * (μ_r - μ0_r) -
            tr(l0_r * Σ_r)
        ) # TODO: in Matlab implementation this is wrong because dim_r should be adjusted after pruning
    log_p_prec =
        a0 * log(b0) - loggamma(a0) + (a0 - 1.0) * (digamma(a_r) - log(b_r)) - b0 * τ_r
    log_p_z = sum(
        log.(1.0 .- p0[z_idx]) .+ z_r[z_idx] .* log.(p0[z_idx] ./ (1.0 .- p0[z_idx]))
    )
    log_q_weight = 0.5 * (logdet(Σ_r) + dim_r * (1.0 + log(2π)))
    log_q_prec = a_r - log(b_r) + loggamma(a_r) + (1.0 - a_r) * digamma(a_r)
    log_q_z = sum(
        -(1.0 .- z_r[z_idx]) .* log.(1.0 .- z_r[z_idx]) .- z_r[z_idx] .* log.(z_r[z_idx])
    )

    return log_lik +
           log_p_weight +
           log_p_prec +
           log_p_z +
           log_q_weight +
           log_q_prec +
           log_q_z
end
