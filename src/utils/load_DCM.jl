"""
    load_DCM(path::String;verbose=true)

Load a DCM from a specified path and return either a linear, bilinear or nonlinear DCM
based on the contents of the file.

# Arguments
- `path::String`: Path to the file to load. Can end with .mat for Matlab files or .jls for data serialized by Julia.
- `verbose::Bool`: Verbosity
- `dcm_key::String`: It is assumed that the variable in MATLAB that was saved as a DCM.mat file was called "DCM". If this is not the case, set `dcm_key` to the variable name that was used (applies only when loading .mat files).

# Output
- `dcm`: DCM struct (can be a linear, bilinear on nonlinear DCM)

# Examples
```julia-repl
julia> dcm = load_DCM("myDCM.jls";verbose=false)
```
"""
function load_DCM(path::String; verbose=true, dcm_key="DCM")
    extension = splitext(path)[end]
    if extension == ".mat"
        return load_MAT_DCM(path; verbose=verbose, dcm_key=dcm_key)
    elseif extension == ".jls"
        return deserialize(path)
    else
        error("Invalid file extension.")
    end
end

"""
$(SIGNATURES)

Save a DCM struct either as .mat or .jls file.

# Arguments
- `path::String`: Path defining where to save file. If ends with .mat -> saves as Matlab file.
    If ends with .jls -> Serializes using Julias Serialization.
- `dcm<:DCM`: DCM struct (can be a linear, bilinear on nonlinear DCM)
# Examples
```julia-repl
julia> save_DCM("./myDCM.mat",dcm)
```
"""
function save_DCM(path::String, dcm::T) where {T<:DCM}
    extension = splitext(path)[end]
    if extension == ".mat"
        save_MAT_DCM(path, dcm)
    elseif extension == ".jls"
        serialize(path, dcm)
    else
        error("Invalid file extension.")
    end
end

"""
    load_MAT_DCM(path;dcm)

Loads a DCM struct in .mat format assuming the SPM naming convention of fields.

# Arguments
- `path::String`: Path of the file to load. Needs to end with .mat
- `verbose::Bool`: Verbosity

# Output
- `dcm::DCM`: DCM struct. Detects automatically if linear, bi-linear or non-linear DCM
based on dcm.b and dcm.d field.
"""
function load_MAT_DCM(path::String; verbose=true, dcm_key="DCM")
    file = matopen(path)
    DCM_mat = read(file, dcm_key)
    close(file)

    # try to extract input information
    if haskey(DCM_mat, "U")
        try
            if haskey(DCM_mat["U"], "name")
                U = InputU(DCM_mat["U"]["u"], DCM_mat["U"]["dt"], vec(DCM_mat["U"]["name"]))
            else
                U = InputU(DCM_mat["U"]["u"], DCM_mat["U"]["dt"])
            end
        catch e
            if isa(e, KeyError)
                error("U struct does not have necessary fields.")
            else
                rethrow(e)
            end
        end
    else
        U = nothing
        if verbose
            @info "No input found-> constructing resting-state DCM."
        end
    end

    # try to extract BOLD signal
    if haskey(DCM_mat, "Y")
        try
            if all(size(DCM_mat["Y"]["y"]) .== 0)
                Y = BoldY(nothing, DCM_mat["Y"]["dt"], nothing)
            else
                if haskey(DCM_mat["Y"], "name")
                    Y = BoldY(
                        DCM_mat["Y"]["y"], DCM_mat["Y"]["dt"], vec(DCM_mat["Y"]["name"])
                    )
                else
                    Y = BoldY(DCM_mat["Y"]["y"], DCM_mat["Y"]["dt"])
                end
            end
        catch e
            if isa(e, KeyError)
                error("Y struct does not have necessary fields.")
            else
                rethrow(e)
            end
        end
    else
        Y = nothing
    end

    # try to extract confounds
    if haskey(DCM_mat, "Conf")
        try
            if haskey(DCM_mat["Conf"], "name")
                Conf = Confound(DCM_mat["Conf"]["X0"], vec(DCM_mat["Conf"]["name"]))
            else
                Conf = Confound(DCM_mat["Conf"]["X0"])
            end
        catch e
            if isa(e, KeyError)
                #@info "Conf struct does not have any fields."
                Conf = nothing
            else
                rethrow(e)
            end
        end
    else
        Conf = nothing
    end

    #---------------------------------------------------------------------------------------
    # check if non-linear DCM
    #---------------------------------------------------------------------------------------
    if haskey(DCM_mat, "d")
        if sum(DCM_mat["d"]) != 0 && all(size(DCM_mat["d"]) .!= 0)
            if !haskey(DCM_mat, "Ep")
                if verbose
                    @info "Ep field missing, using default values."
                end
                Tp = TrueParamNonLinear(
                    BitMatrix(DCM_mat["a"]),
                    BitArray(DCM_mat["b"]),
                    BitMatrix(DCM_mat["c"]),
                    BitArray(DCM_mat["d"]);
                    sample=false,
                )
            else
                if DCM_mat["Ep"]["transit"] isa Vector
                    Tp = TrueParamNonLinear(
                        DCM_mat["Ep"]["A"],
                        DCM_mat["Ep"]["B"],
                        DCM_mat["Ep"]["C"],
                        DCM_mat["Ep"]["D"],
                        DCM_mat["Ep"]["transit"],
                        DCM_mat["Ep"]["decay"],
                        only(DCM_mat["Ep"]["epsilon"]),
                    )
                else
                    Tp = TrueParamNonLinear(
                        DCM_mat["Ep"]["A"],
                        DCM_mat["Ep"]["B"],
                        DCM_mat["Ep"]["C"],
                        DCM_mat["Ep"]["D"],
                        vec(Matrix(DCM_mat["Ep"]["transit"])),
                        vec(Matrix(DCM_mat["Ep"]["decay"])),
                        only(DCM_mat["Ep"]["epsilon"]),
                    )
                end
            end

            if verbose
                @info "Found non-linear DCM."
            end
            return NonLinearDCM(
                DCM_mat["a"],
                DCM_mat["b"],
                DCM_mat["c"],
                DCM_mat["d"],
                Int64(DCM_mat["v"]),
                Int64(DCM_mat["n"]),
                U,
                Y,
                Tp,
                Conf,
            )
        end
    end

    #---------------------------------------------------------------------------------------
    # check if bi-linear DCM
    #---------------------------------------------------------------------------------------
    if haskey(DCM_mat, "b") && sum(DCM_mat["b"]) != 0.0 && all(size(DCM_mat["b"]) .!= 0)
        if !haskey(DCM_mat, "Ep")
            if verbose
                @info "Ep field missing, using default values."
            end
            Tp = TrueParamBiLinear(
                BitMatrix(DCM_mat["a"]),
                BitArray(DCM_mat["b"]),
                BitMatrix(DCM_mat["c"]);
                sample=false,
            )
        else
            if DCM_mat["Ep"]["transit"] isa Vector
                Tp = TrueParamBiLinear(
                    DCM_mat["Ep"]["A"],
                    DCM_mat["Ep"]["B"],
                    DCM_mat["Ep"]["C"],
                    DCM_mat["Ep"]["transit"],
                    DCM_mat["Ep"]["decay"],
                    only(DCM_mat["Ep"]["epsilon"]),
                )
            else
                Tp = TrueParamBiLinear(
                    DCM_mat["Ep"]["A"],
                    DCM_mat["Ep"]["B"],
                    DCM_mat["Ep"]["C"],
                    vec(Matrix(DCM_mat["Ep"]["transit"])),
                    vec(Matrix(DCM_mat["Ep"]["decay"])),
                    only(DCM_mat["Ep"]["epsilon"]),
                )
            end
        end

        if verbose
            @info "Found bi-linear DCM."
        end
        return BiLinearDCM(
            DCM_mat["a"],
            DCM_mat["b"],
            DCM_mat["c"],
            Int64(DCM_mat["v"]),
            Int64(DCM_mat["n"]),
            U,
            Y,
            Tp,
            Conf,
        )
    end

    #---------------------------------------------------------------------------------------
    # if non of the above -> must be a linear DCM
    #---------------------------------------------------------------------------------------
    if !haskey(DCM_mat, "Ep")
        if verbose
            @info "Ep field missing, using default values."
        end
        Tp = TrueParamLinear(BitMatrix(DCM_mat["a"]), BitMatrix(DCM_mat["c"]); sample=false)
    else
        try
            if DCM_mat["Ep"]["transit"] isa Vector
                Tp = TrueParamLinear(
                    DCM_mat["Ep"]["A"],
                    DCM_mat["Ep"]["C"],
                    DCM_mat["Ep"]["transit"],
                    DCM_mat["Ep"]["decay"],
                    only(DCM_mat["Ep"]["epsilon"]),
                )
            else
                Tp = TrueParamLinear(
                    DCM_mat["Ep"]["A"],
                    DCM_mat["Ep"]["C"],
                    vec(Matrix(DCM_mat["Ep"]["transit"])),
                    vec(Matrix(DCM_mat["Ep"]["decay"])),
                    only(DCM_mat["Ep"]["epsilon"]),
                )
            end
        catch e
            if isa(e, KeyError)
                error("Ep struct does not have necessary fields.")
            else
                rethrow(e)
            end
        end
    end

    if verbose
        @info "Found linear DCM."
    end
    return LinearDCM(
        DCM_mat["a"], DCM_mat["c"], Int64(DCM_mat["v"]), Int64(DCM_mat["n"]), U, Y, Tp, Conf
    )
end

"""
    save_MAT_DCM(path,dcm)

Saves a DCM struct in .mat format using the SPM naming convention.

# Arguments
- `path::String`: Path specifying where to save file. Needs to end with .mat
- `dcm::DCM`: DCM struct to be saved (can be linear, bi-linear or non-linear DCM)
"""
function save_MAT_DCM(path::String, dcm::T) where {T<:DCM}
    dcm_dict = construct_dict(dcm)
    return matwrite(path, Dict("DCM" => dcm_dict))
end

function getfield_special(obj::T, key::Symbol) where {T<:DCM}
    prop = getproperty(obj, key)
    if prop isa BoldY ||
        prop isa InputU ||
        prop isa TrueParamLinear ||
        prop isa TrueParamBiLinear ||
        prop isa TrueParamNonLinear ||
        prop isa Confound
        return Dict(if getfield(prop, k) isa Nothing
            String(k) => empty_field(k)
        else
            String(k) => getfield(prop, k)
        end for k in propertynames(prop))
    elseif prop isa Nothing
        return Dict()
    else
        return getfield(obj, key)
    end
end

function empty_field(k::Symbol)
    if k === :y
        return Array{Float64}(undef, 0, 0)
    elseif k === :name
        return Array{String}(undef, 0)
    end
end

"""
$(SIGNATURES)

Loads an example DCM (linear DCM with 50 regions). The network architecture is based on the
S50 structure introduced in [Smith et al. (2011). Network modelling methods for FMRI. NeuroImage.](https://pubmed.ncbi.nlm.nih.gov/20817103/)
This DCM can be used to generate synthetic data.

# Examples
```julia-repl
julia> dcm = load_example_DCM()
julia> y_noise, _, _, _ = generate_BOLD(dcm;SNR=10)
```
"""
function load_example_DCM()
    path = joinpath(artifact"rDCMdata", "DCM_LargeScaleSmith_model1.mat")
    return load_DCM(path; verbose=false)
end

function load_DCM_bilinear()
    path = joinpath(artifact"rDCMdata", "DCM_bilinear.mat")
    return load_DCM(path; verbose=false)
end

function load_DCM_nonlinear()
    path = joinpath(artifact"rDCMdata", "DCM_nonlinear.mat")
    return load_DCM(path; verbose=false)
end
