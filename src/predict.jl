"""
$(TYPEDSIGNATURES)

Predict BOLD signal based on parameter estimates after model inversion.

# Arguments
- `rdcm`: rDCM struct that was used to estimate parameters.
- `output`: Output of model inversion.

# Output
- `y_pred::Matrix{Float64}`: predicted BOLD signal
"""
function predict(rdcm::T1, output::T2) where {T1<:RDCM,T2<:ModelOutput}
    out_dcm = LinearDCM(rdcm, output)
    # SNR value has no meaning because we save only noise free prediction
    _, y_pred, _, _ = generate_BOLD(out_dcm; SNR=1)
    return y_pred
end
