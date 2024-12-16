
function get_prior_stats(a::BitMatrix, c::BitMatrix)
    nr = size(a, 1)

    fac = 8 # for endogenous DCMs it would be 128, TODO: check if we need that case

    a_new = a - diagm(diag(a))

    # prior expectation
    #pE_A = a_new./(64nr) - I/2 # original SPM variant
    pE_A = zeros(nr, nr) - I / 2 # rDCM variant
    pE_C = zeros(Float64, size(c))

    # prior mean
    pE_transit = zeros(Float64, nr)
    pE_decay = zeros(Float64, nr)
    pE_epsilon = 0.0

    # prior covariance
    pC_A = (a_new .* fac) ./ nr + I ./ (8nr)
    pC_C = zeros(Float64, size(c))
    pC_C[c] .= 1.0

    pC_transit = zeros(Float64, nr) .+ exp(-6)
    pC_decay = zeros(Float64, nr) .+ exp(-6)
    pC_epsilon = exp(-6)

    pE = priorMeanLinear(pE_A, pE_C, pE_transit, pE_decay, pE_epsilon)
    pC = priorCovLinear(pC_A, pC_C, pC_transit, pC_decay, pC_epsilon)

    return PriorDCMLinear(pE, pC)
end

function get_prior_stats(a::BitMatrix, b::BitArray{3}, c::BitMatrix)
    nr = size(a, 1)

    fac = 8 # for endogenous DCMs it would be 128, TODO: check if we need that case

    a_new = a - diagm(diag(a))

    # prior expectation
    #pE_A = a_new./(64nr) - I/2 # original SPM variant
    pE_A = zeros(nr, nr) - I / 2 # rDCM variant
    pE_B = zeros(Float64, size(b))
    pE_C = zeros(Float64, size(c))

    # prior mean
    pE_transit = zeros(Float64, nr)
    pE_decay = zeros(Float64, nr)
    pE_epsilon = 0.0

    # prior covariance
    pC_A = (a_new .* fac) ./ nr + I ./ (8nr)
    pC_B = zeros(Float64, size(b))
    pC_C = zeros(Float64, size(c))
    pC_B[b] .= 1.0
    pC_C[c] .= 1.0

    pC_transit = zeros(Float64, nr) .+ exp(-6)
    pC_decay = zeros(Float64, nr) .+ exp(-6)
    pC_epsilon = exp(-6)

    pE = priorMeanBiLinear(pE_A, pE_B, pE_C, pE_transit, pE_decay, pE_epsilon)
    pC = priorCovBiLinear(pC_A, pC_B, pC_C, pC_transit, pC_decay, pC_epsilon)

    return PriorDCMBiLinear(pE, pC)
end

function get_prior_stats(a::BitMatrix, b::BitArray{3}, c::BitMatrix, d::BitArray{3})
    nr = size(a, 1)

    fac = 8 # for endogenous DCMs it would be 128, TODO: check if we need that case

    a_new = a - diagm(diag(a))

    # prior expectation
    #pE_A = a_new./(64nr) - I/2 # original SPM variant
    pE_A = zeros(nr, nr) - I / 2 # rDCM variant
    pE_B = zeros(Float64, size(b))
    pE_C = zeros(Float64, size(c))
    pE_D = zeros(Float64, size(d))

    # prior mean
    pE_transit = zeros(Float64, nr)
    pE_decay = zeros(Float64, nr)
    pE_epsilon = 0.0

    # prior covariance
    pC_A = (a_new .* fac) ./ nr + I ./ (8nr)
    pC_B = zeros(Float64, size(b))
    pC_C = zeros(Float64, size(c))
    pC_D = zeros(Float64, size(d))
    pC_B[b] .= 1.0
    pC_C[c] .= 1.0
    pC_D[d] .= 1.0

    pC_transit = zeros(Float64, nr) .+ exp(-6)
    pC_decay = zeros(Float64, nr) .+ exp(-6)
    pC_epsilon = exp(-6)

    pE = priorMeanNonLinear(pE_A, pE_B, pE_C, pE_D, pE_transit, pE_decay, pE_epsilon)
    pC = priorCovNonLinear(pC_A, pC_B, pC_C, pC_D, pC_transit, pC_decay, pC_epsilon)

    return PriorDCMNonLinear(pE, pC)
end

function sample_from_prior(a::BitMatrix, c::BitMatrix; fixHRF=true, rng=Xoshiro())
    prior = rDCM.get_prior_stats(a, c)
    return sample_from_prior(a, c, prior; fixHRF=fixHRF, rng=rng)
end

function sample_from_prior(
    a::BitMatrix, c::BitMatrix, prior::PriorDCMLinear; fixHRF=true, rng=Xoshiro()
)

    # prepare mean and variance for connectivity parameters
    mu0 = [prior.pE.A[a]; prior.pE.C[c]]
    sigma0 = [prior.pC.A[a]; prior.pC.C[c]]
    # sigma0 = diagm([prior.pC.A[a]; prior.pC.C[c]])

    A_endIdx = sum(a)

    for _ in 1:50
        # mu = rand(rng, MvNormal(mu0, sigma0))
        mu = similar(mu0)
        for i in eachindex(mu)
            mu[i] = rand(rng, Normal(mu0[i], sigma0[i]))
        end

        A = zeros(size(a))
        A[a] = mu[1:A_endIdx]

        # check if system is stable
        if maximum(real(eigvals(A))) < 0
            C = zeros(size(c))
            C[c] = mu[(A_endIdx + 1):end]

            if fixHRF
                transit = prior.pE.transit
                decay = prior.pE.decay
                epsilon = prior.pE.epsilon
            else
                # sample HRF parameters
                transit = rand(rng, MvNormal(prior.pE.transit, diagm(prior.pC.transit)))
                decay = rand(rng, MvNormal(prior.pE.decay, diagm(prior.pC.decay)))
                epsilon = rand(rng, Normal(prior.pE.epsilon, prior.pC.epsilon))
            end

            return A, C, transit, decay, epsilon
        end
    end
    return error("Not able so sample values such that system is stable.")
end

function sample_from_prior(
    a::BitMatrix,
    b::BitArray{3},
    c::BitMatrix,
    prior::PriorDCMBiLinear;
    fixHRF=true,
    rng=Xoshiro(),
)

    # prepare mean and variance for connectivity parameters
    mu0 = [prior.pE.A[a]; prior.pE.B[b]; prior.pE.C[c]]
    sigma0 = diagm([prior.pC.A[a]; prior.pC.B[b]; prior.pC.C[c]])

    #rng = Xoshiro(seed)

    A_endIdx = sum(a)
    len_b = sum(b)

    for _ in 1:50
        mu = rand(rng, MvNormal(mu0, sigma0))

        A = zeros(size(a))
        A[a] = mu[1:A_endIdx]

        # check if system is stable
        if maximum(real(eigvals(A))) < 0
            B = zeros(size(b))
            C = zeros(size(c))
            B[b] = mu[(A_endIdx + 1):(A_endIdx + len_b)]
            C[c] = mu[(A_endIdx + 1 + len_b):end]

            if fixHRF
                transit = prior.pE.transit
                decay = prior.pE.decay
                epsilon = prior.pE.epsilon
            else
                # sample HRF parameters
                transit = rand(rng, MvNormal(prior.pE.transit, diagm(prior.pC.transit)))
                decay = rand(rng, MvNormal(prior.pE.decay, diagm(prior.pC.decay)))
                epsilon = rand(rng, Normal(prior.pE.epsilon, prior.pC.epsilon))
            end

            return A, B, C, transit, decay, epsilon
        end
    end
    return error("Not able so sample values such that system is stable.")
end

function sample_from_prior(
    a::BitMatrix,
    b::BitArray{3},
    c::BitMatrix,
    d::BitArray{3},
    prior::PriorDCMNonLinear;
    fixHRF=true,
    rng=Xoshiro(),
)

    # prepare mean and variance for connectivity parameters
    mu0 = [prior.pE.A[a]; prior.pE.B[b]; prior.pE.C[c]; prior.pE.D[d]]
    sigma0 = diagm([prior.pC.A[a]; prior.pC.B[b]; prior.pC.C[c]; prior.pC.D[d]])

    #rng = Xoshiro(seed)

    A_endIdx = sum(a)
    len_b = sum(b)
    len_c = sum(c)

    for _ in 1:50
        mu = rand(rng, MvNormal(mu0, sigma0))

        A = zeros(size(a))
        A[a] = mu[1:A_endIdx]

        # check if system is stable
        if maximum(real(eigvals(A))) < 0
            B = zeros(size(b))
            C = zeros(size(c))
            D = zeros(size(d))
            B[b] = mu[(A_endIdx + 1):(A_endIdx + len_b)]
            C[c] = mu[(A_endIdx + 1 + len_b):(A_endIdx + len_b + len_c)]
            D[d] = mu[(A_endIdx + 1 + len_b + len_c):end]

            if fixHRF
                transit = prior.pE.transit
                decay = prior.pE.decay
                epsilon = prior.pE.epsilon
            else
                # sample HRF parameters
                transit = rand(rng, MvNormal(prior.pE.transit, diagm(prior.pC.transit)))
                decay = rand(rng, MvNormal(prior.pE.decay, diagm(prior.pC.decay)))
                epsilon = rand(rng, Normal(prior.pE.epsilon, prior.pC.epsilon))
            end

            return A, B, C, D, transit, decay, epsilon
        end
    end
    return error("Not able so sample values such that system is stable.")
end

function get_priors(rdcm::RigidRdcm, nu::Int64)
    prior = get_prior_stats(rdcm.a, rdcm.c)

    # prior mean
    m0 = [prior.pE.A prior.pE.C]

    # prior precision
    pC_A = 1 ./ prior.pC.A
    pC_C = 1 ./ prior.pC.C

    # prior precision of baseline
    #if sum(pC_C[:, end]) ≠ 0
    #    pC_C[:, end] .= 1.0e-8
    #end

    # prior precision
    l0 = [pC_A pC_C]

    # setting priors on noise precision
    a0 = 2.0
    b0 = 1.0

    return m0, l0, a0, b0
end

function get_priors(rdcm::BiLinearRigidRdcm, nu::Int64)
    prior = get_prior_stats(rdcm.a, rdcm.b, rdcm.c)

    nr = size(prior.pE.A,1)

    # prior mean
    m0 = [prior.pE.A reshape(prior.pE.B,nr,nr*nu) prior.pE.C]

    # prior precision
    pC_A = 1 ./ prior.pC.A
    pC_B = 1 ./ prior.pC.B
    pC_C = 1 ./ prior.pC.C

    # prior precision of baseline
    #if sum(pC_C[:, end]) ≠ 0
    #    pC_C[:, end] .= 1.0e-8
    #end

    # prior precision
    l0 = [pC_A reshape(pC_B,nr,nr*nu) pC_C]

    # setting priors on noise precision
    a0 = 2.0
    b0 = 1.0

    return m0, l0, a0, b0
end

function get_priors(rdcm::T, conf_weights::Union{BitVector,BitMatrix}) where {T<:RDCM}
    dcm = copy(rdcm)
    nu = size(dcm.c,2)
    dcm.c = [dcm.c conf_weights]
    return get_priors(dcm, nu)
end

function get_priors(rdcm::SparseRdcm, nu::Int64)
    prior = get_prior_stats(BitMatrix(ones(size(rdcm.a))), BitMatrix(ones(size(rdcm.c))))

    A = prior.pE.A
    # set the prior mean of endogenous parameters to zero
    A = zeros(size(A)) + diagm(diag(A))

    # prior mean
    m0 = [A prior.pE.C]

    # prior precision
    pC_A = 1 ./ prior.pC.A
    pC_C = 1 ./ prior.pC.C

    # TODO: ask Stefan why not set prior precision for baseline like for rigid rdcm

    # prior precision
    l0 = [pC_A pC_C]

    # setting priors on noise precision
    a0 = 2.0
    b0 = 1.0

    return m0, l0, a0, b0
end
