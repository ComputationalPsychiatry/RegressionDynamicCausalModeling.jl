"""
    RegressionDynamicCausalModeling

A Julia package for estimating whole-brain effective connectivity using the regression
dynamic causal modeling (rDCM) framework.

The alias `rDCM` is exported for the package name.
"""
module RegressionDynamicCausalModeling

const rDCM = RegressionDynamicCausalModeling

using Distributions
using DocStringExtensions
using FFTW: rfft, irfft, fft, ifft
using LinearAlgebra: tr, inv, logdet, diagm, diag, I, diagind, eigvals, isposdef, Hermitian
using MAT: matopen, matwrite
using LazyArtifacts
using PrecompileTools: @compile_workload, @setup_workload
using Printf
using ProgressMeter
using Random
using Scratch
using Serialization: serialize, deserialize
using SparseArrays: spzeros, SparseMatrixCSC, SparseArrays
using SpecialFunctions: digamma, loggamma
using Statistics

const euler_integration_bin = @get_scratch!("euler_int_bin") # location of binary for Euler integration
const tmpdir = @get_scratch!("tmp")

"""

$(TYPEDEF)

Many models and functions accept a `verbosity` parameter.

Choose between: `NONE`, `STD` [default], `FULL`.
"""
@enum Verbosity NONE = 0 STD = 1 FULL = 2

"""
    const FIXEDSEED

Fixed seed to allow reproducible results.
"""
const FIXEDSEED = 42

include("structs.jl")
include("constructors.jl")
include("create_regressors.jl")
include("utils/load_DCM.jl")
include("utils/dcm_euler_integration.jl")
include("utils/dcm_print.jl")
include("utils/SPM_compat.jl")
include("generate_BOLD.jl")
include("get_priors.jl")
include("rigid_inversion.jl")
include("invert_rDCM.jl")
include("predict.jl")
include("sparse_inversion.jl")
include("precompile.jl")

export rDCM

# IO
export load_DCM, save_DCM, load_example_DCM, export_to_SPM
# Input, BOLD and Confound
export InputU, BoldY, Confound
# DCM structs
export DCM, LinearDCM, BiLinearDCM, NonLinearDCM, RDCM, RigidRdcm, SparseRdcm
# inversion
export Options, RigidInversionParams, SparseInversionParams, invert
# output
export RigidOutput, SparseOutput
# prediction
export predict
# misc
export generate_BOLD, sample_from_prior

end
