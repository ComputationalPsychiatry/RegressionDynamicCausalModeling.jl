naming_rDCM_to_SPM = Base.ImmutableDict("scans" => "v", "nr" => "n")

"""
$(SIGNATURES)

Export a DCM with output after model inversion as an SPM compatible .mat file.

# Arguments
- `path::String`: Path defining where to save file. Needs to end with .mat
- `rdcm::RigidRdcm`: An rDCM model.
- `output::RigidOutput`: Output after model inversion.
# Examples
```julia-repl
julia> save_DCM("DCM.mat",rdcm,output)
```

!!! info
    See [SPM compatibility](examples/SPM_compat.md) for limitations of this functionality.
"""
function export_to_SPM(path::String, rdcm::RigidRdcm, output::RigidOutput)
    dcm = LinearDCM(rdcm, output)
    dcm_dict = construct_dict(dcm)

    # prior mean and covariance
    m0, l0, _, _ = get_priors(rdcm)
    spm_pE = Dict(
        "A" => m0[1:size(dcm.a, 1), 1:size(dcm.a, 2)],
        "C" => m0[1:size(dcm.c, 1), (size(dcm.a, 2) + 1):end],#TODO check this
    )

    pc_A = 1 ./ l0[1:size(dcm.a, 1), 1:size(dcm.a, 2)]
    pc_C = 1 ./ l0[1:size(dcm.c, 1), (size(dcm.a, 2) + 1):end] # TODO check this

    n_param = size(dcm.a, 1)*size(dcm.a, 2) + size(dcm.c, 1)*size(dcm.c, 2)
    spm_pC = spzeros(n_param, n_param)
    spm_pC[diagind(spm_pC)] .= [pc_A[:]; pc_C[:]]
    # model settings
    spm_M = Dict(
        "IS" => output.inversion, "pE" => spm_pE, "pC" => spm_pC
    )

    # error covariance
    Ce = 1 ./ (output.a_all ./ output.b_all)
    merge!(dcm_dict, Dict("Ce" => Ce))

    merge!(dcm_dict, Dict("M" => spm_M))
    merge!(dcm_dict, Dict("F" => output.F)) #save also negative free energy
    return matwrite(path, Dict("DCM" => dcm_dict))
end

function construct_dict(dcm::T) where {T<:DCM}
    dcm_dict = Dict{String,Any}(
        String(SPM_mapping(key)) => getfield_special(dcm, key) for key in propertynames(dcm)
    )
    if dcm isa NonLinearDCM
        spm_options = Dict(
            "nonlinear" => 1,
            "two_state" => 0,
            "stochastic" => 0,
            "centre" => 0,
            "induced" => 0,
        )
        merge!(dcm_dict, Dict("options" => spm_options))
    else
        spm_options = Dict(
            "nonlinear" => 0,
            "two_state" => 0,
            "stochastic" => 0,
            "centre" => 0,
            "induced" => 0,
        )
        merge!(dcm_dict, Dict("options" => spm_options))
        merge!(dcm_dict, Dict("d" => Array{Float64}(undef, dcm.nr, dcm.nr, 0)))
        if dcm isa LinearDCM
            merge!(dcm_dict, Dict("b" => zeros(dcm.nr, dcm.nr, size(dcm.c, 2))))
        end
    end
    return dcm_dict
end

function SPM_mapping(key::Symbol)
    if haskey(naming_rDCM_to_SPM, String(key))
        return get(naming_rDCM_to_SPM, String(key), Inf)
    else
        return key
    end
end
