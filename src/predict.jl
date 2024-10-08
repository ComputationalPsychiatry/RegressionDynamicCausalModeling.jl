"""
$(TYPEDSIGNATURES)

Predict BOLD signal based on parameter estimates after model inversion.

# Arguments
- `rdcm`: rDCM struct that was used to estimate parameters.
- `output`: Output of model inversion.

# Output
- `y_pred::Matrix{Float64}`: predicted BOLD signal

!!! warning
    The predicted signal is generated by integrating a linear DCM based on the estimated
    parameters. Therefore, it is not possible to generate predictions for resting-state fMRI
    since there is no driving input that would elicit activity.

"""
function predict(rdcm::T1, output::T2) where {T1<:RDCM,T2<:ModelOutput}
    out_dcm = LinearDCM(rdcm, output)
    if sum(out_dcm.c) == 0 || isnothing(out_dcm.U) || size(out_dcm.U.u, 1) == 0
        error("Cannot generate data from resting-state DCM.")
    end
    # SNR value has no meaning because we only save noise free prediction
    _, y_pred, _, _ = generate_BOLD(out_dcm; SNR=1)
    return y_pred
end
