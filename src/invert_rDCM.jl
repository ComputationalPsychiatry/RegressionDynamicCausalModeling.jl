"""
$(TYPEDSIGNATURES)

Invert an rDCM model specified by `rdcm` with fixed network architecture using specific inversion settings defined in `opt`.
"""
function invert(rdcm::RigidRdcm, opt::Options)
    if opt.verbose > 0
        println("Running model inversion (rigid rDCM)") # COV_EXCL_LINE
    end

    X, Y = create_regressors!(rdcm, opt.rng) # The mathematically correct version (where we don't have duplicated data)

    return rigid_inversion(rdcm, X, Y, opt) # new mathematically correct version
end

"""
$(TYPEDSIGNATURES)

Invert a sparse rDCM model specified by `rdcm` using specific inversion settings defined in `opt`.
"""
function invert(rdcm::SparseRdcm, opt::Options)
    if opt.verbose > 0
        println("Running model inversion (sparse rDCM)") # COV_EXCL_LINE
    end

    X, Y = create_regressors!(rdcm, opt.rng) # The mathematically correct version (where we don't have duplicated data)

    return sparse_inversion(rdcm, X, Y, opt) # new mathematically correct version
end

function invert(rdcm::BiLinearRigidRdcm, opt::Options)
    if opt.verbose > 0
        println("Running model inversion (rigid rDCM)") # COV_EXCL_LINE
    end

    X, Y = create_regressors!(rdcm, opt.rng)

    return rigid_inversion(rdcm, X, Y, opt)
end
