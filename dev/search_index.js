var documenterSearchIndex = {"docs":
[{"location":"api/#API-reference","page":"API reference","title":"API reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"RegressionDynamicCausalModeling","category":"page"},{"location":"api/#RegressionDynamicCausalModeling","page":"API reference","title":"RegressionDynamicCausalModeling","text":"RegressionDynamicCausalModeling\n\nA Julia package for estimating whole-brain effective connectivity using the regression dynamic causal modeling (rDCM) framework.\n\nThe alias rDCM is exported for the package name.\n\n\n\n\n\n","category":"module"},{"location":"api/#DCMs","page":"API reference","title":"DCMs","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"DCM\nLinearDCM\nLinearDCM(a,c,scans,nr,U,Y,Ep,Conf)\nLinearDCM(::Matrix{Number},::Matrix{Number},::Int64,::Int64,::InputU,::Union{BoldY,Nothing},::rDCM.TrueParamLinear)\nLinearDCM(::DCM)\nLinearDCM(::RDCM,::rDCM.ModelOutput)\nBiLinearDCM\nBiLinearDCM(a,b,c,scans,nr,U,Y,Ep,Conf)\nNonLinearDCM\nNonLinearDCM(a,b,c,d,scans,nr,U,Y,Ep,Conf)\nInputU\nBoldY","category":"page"},{"location":"api/#RegressionDynamicCausalModeling.DCM","page":"API reference","title":"RegressionDynamicCausalModeling.DCM","text":"DCM\n\nAbstract supertype for a DCM.\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.LinearDCM","page":"API reference","title":"RegressionDynamicCausalModeling.LinearDCM","text":"mutable struct LinearDCM <: DCM\n\nRepresentation of a linear DCM.\n\nFields\n\na::BitMatrix: Binary indicator matrix for endogenous connectivity\nc::BitMatrix: Binary indicator matrix for driving inputs\nscans::Int64: number of data points per region\nnr::Int64\nU::Union{Nothing, InputU}: input structure with information about driving input\nY::Union{Nothing, BoldY}: data structure containing BOLD signal\nEp::RegressionDynamicCausalModeling.TrueParamLinear: connectivity parameters (A and C matrix)\nConf::Union{Nothing, RegressionDynamicCausalModeling.Confound}: confound structure\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.LinearDCM-NTuple{8, Any}","page":"API reference","title":"RegressionDynamicCausalModeling.LinearDCM","text":"LinearDCM(a, c, scans, nr, U, Y, Ep, Conf)\n\n\nConstructor with sanity checks for LinearDCM. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.LinearDCM-Tuple{Matrix{Number}, Matrix{Number}, Int64, Int64, InputU, Union{Nothing, BoldY}, RegressionDynamicCausalModeling.TrueParamLinear}","page":"API reference","title":"RegressionDynamicCausalModeling.LinearDCM","text":"LinearDCM(\n    a::Array{T1<:Number, 2},\n    c::Array{T2<:Number, 2},\n    scans::Int64,\n    nr::Int64,\n    U::InputU,\n    Y::Union{Nothing, BoldY},\n    Ep::RegressionDynamicCausalModeling.TrueParamLinear\n) -> LinearDCM\n\n\nCreate linear DCM.\n\nArguments\n\na::Matrix: Binary indicator matrix for endogenous connectivity\nc::Matrix: Binary indicator matrix for driving inputs\nscans::Int64: Number of data points per region\nnr::Int64: Number of regions\nU::InputU: Input structure\nY::Union{BoldY,Nothing}: Data structure containing BOLD signal, can be nothing\nEP::TrueParamLinear: Connectivity parameters containing A and C matrix\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.LinearDCM-Tuple{DCM}","page":"API reference","title":"RegressionDynamicCausalModeling.LinearDCM","text":"LinearDCM(dcm::DCM) -> LinearDCM\n\n\nCreate linear DCM from a bilinear or nonlinear DCM dcm.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.LinearDCM-Tuple{RDCM, RegressionDynamicCausalModeling.ModelOutput}","page":"API reference","title":"RegressionDynamicCausalModeling.LinearDCM","text":"LinearDCM(\n    rdcm::RDCM,\n    output::RegressionDynamicCausalModeling.ModelOutput\n) -> LinearDCM\n\n\nCreate linear DCM from an rDCM structure rdcm and output from model inversion output.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.BiLinearDCM","page":"API reference","title":"RegressionDynamicCausalModeling.BiLinearDCM","text":"mutable struct BiLinearDCM <: DCM\n\nRepresenetation of a bi-linear DCM.\n\nFields\n\na::BitMatrix: Binary indicator matrix for endogenous connectivity\nb::BitArray{3}: Binary indicator matrix for bi-linear dynamics.\nc::BitMatrix: Binary indicator matrix for driving inputs\nscans::Int64: number of data points per region\nnr::Int64\nU::Union{Nothing, InputU}: input structure with information about driving input\nY::Union{Nothing, BoldY}: data structure containing BOLD signal\nEp::RegressionDynamicCausalModeling.TrueParamBiLinear: connectivity parameters (A, B and C matrix)\nConf::Union{Nothing, RegressionDynamicCausalModeling.Confound}: confound structure\n\nwarning: Warning\nWhile this package allows to simulate data for bi-linear DCMs, the current rDCM implementation can only invert linear DCMs.\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.BiLinearDCM-NTuple{9, Any}","page":"API reference","title":"RegressionDynamicCausalModeling.BiLinearDCM","text":"BiLinearDCM(a, b, c, scans, nr, U, Y, Ep, Conf)\n\n\nConstructor with sanity checks for BiLinearDCM. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.NonLinearDCM","page":"API reference","title":"RegressionDynamicCausalModeling.NonLinearDCM","text":"mutable struct NonLinearDCM <: DCM\n\nRepresentation of a non-linear DCM.\n\nFields\n\na::BitMatrix: Binary indicator matrix for endogenous connectivity\nb::BitArray{3}: Binary indicator matrix for bi-linear dynamics.\nc::BitMatrix: Binary indicator matrix for driving inputs\nd::BitArray{3}: Binary indicator matrix for non-linear dynamics\nscans::Int64: number of data points per region\nnr::Int64\nU::Union{Nothing, InputU}: input structure with information about driving input\nY::Union{Nothing, BoldY}: data structure containing BOLD signal\nEp::RegressionDynamicCausalModeling.TrueParamNonLinear: connectivity parameters (A, B, C and D matrix)\nConf::Union{Nothing, RegressionDynamicCausalModeling.Confound}: confound structure\n\nwarning: Warning\nWhile this package allows to simulate data for non-linear DCMs, the current rDCM implementation can only invert linear DCMs.\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.NonLinearDCM-NTuple{10, Any}","page":"API reference","title":"RegressionDynamicCausalModeling.NonLinearDCM","text":"NonLinearDCM(a, b, c, d, scans, nr, U, Y, Ep, Conf)\n\n\nConstructor with sanity checks for NonLinearDCM. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.InputU","page":"API reference","title":"RegressionDynamicCausalModeling.InputU","text":"struct InputU\n\nInput structure for a DCM.\n\nFields\n\nu::Matrix{Float64}: Driving input\ndt::Float64: Sampling interval\nname::Vector{String}: Input names\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.BoldY","page":"API reference","title":"RegressionDynamicCausalModeling.BoldY","text":"mutable struct BoldY\n\nStruct for BOLD signal.\n\nFields\n\ny::Union{Nothing, Matrix{Float64}}: BOLD signal\ndt::Float64: Sampling interval\nname::Union{Nothing, Vector{String}}: Brain region names\n\n\n\n\n\n","category":"type"},{"location":"api/#rDCMs","page":"API reference","title":"rDCMs","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"RDCM\nRigidRdcm\nRigidRdcm(a,c,scans,nr,U,Y,Ep,Conf,hrf)\nRigidRdcm(::LinearDCM)\nRigidOutput\nSparseRdcm\nSparseRdcm(a,c,scans,nr,U,Y,Ep,Conf,hrf,inform_p0,p0)\nSparseRdcm(::LinearDCM;::Bool,::Float64)\nSparseOutput","category":"page"},{"location":"api/#RegressionDynamicCausalModeling.RDCM","page":"API reference","title":"RegressionDynamicCausalModeling.RDCM","text":"RDCM\n\nAbstract supertype for a rigid or sparse rDCM.\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.RigidRdcm","page":"API reference","title":"RegressionDynamicCausalModeling.RigidRdcm","text":"mutable struct RigidRdcm <: RDCM\n\nRepresentation of an rDCM with fixed network architecture.\n\nFields\n\na::BitMatrix: Binary indicator matrix for endogenous connectivity\nc::BitMatrix: Binary indicator matrix for driving inputs\nscans::Int64: Number of data points per region\nnr::Int64\nU::InputU: Input structure with information about driving input\nY::BoldY: Data structure containing BOLD signal\nEp::RegressionDynamicCausalModeling.TrueParamLinear: Connectivity parameters (A and C matrix)\nConf::RegressionDynamicCausalModeling.Confound: Confound structure\nhrf::Vector{Float64}: Hemodynamic response function\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.RigidRdcm-NTuple{9, Any}","page":"API reference","title":"RegressionDynamicCausalModeling.RigidRdcm","text":"RigidRdcm(a, c, scans, nr, U, Y, Ep, Conf, hrf)\n\n\nConstructor with sanity checks for RigidRdcm. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.RigidRdcm-Tuple{LinearDCM}","page":"API reference","title":"RegressionDynamicCausalModeling.RigidRdcm","text":"RigidRdcm(dcm::LinearDCM) -> RigidRdcm\n\n\nConstruct a RigidRdcm based on a linear DCM.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.RigidOutput","page":"API reference","title":"RegressionDynamicCausalModeling.RigidOutput","text":"struct RigidOutput <: RegressionDynamicCausalModeling.ModelOutput\n\nOuput after inversion of RigidRdcm.\n\nFields\n\nF::Float64: Negative free energy of the whole model\nF_r::Vector{Float64}: Region specific negative free energy\niter_all::Vector{Int64}: Number of iterations per region until convergence\na_all::Vector{Float64}: Posterior shape parameter for noise\nb_all::Vector{Float64}: Posterior rate parameter for noise\nm_all::Matrix{Float64}: Posterior mean for connectivity parameters\nΣ_all::Vector{SparseArrays.SparseMatrixCSC{Float64, Int64}}: Posterior covariance for connectivity parameters\ninversion::String: Inversion method\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.SparseRdcm","page":"API reference","title":"RegressionDynamicCausalModeling.SparseRdcm","text":"mutable struct SparseRdcm <: RDCM\n\nRepresentation of a sparse rDCM (sparsity constraints on network architecture).\n\nFields\n\na::BitMatrix: Binary indicator matrix for endogenous connectivity\nc::BitMatrix: Binary indicator matrix for driving inputs\nscans::Int64: Number of data points per region\nnr::Int64\nU::InputU: Input structure with information about driving input\nY::BoldY: Data structure containing BOLD signal\nEp::RegressionDynamicCausalModeling.TrueParamLinear: Connectivity parameters (A and C matrix)\nConf::RegressionDynamicCausalModeling.Confound: Confound structure\nhrf::Vector{Float64}: Hemodynamic response function\ninform_p0::Bool: Inform region specific sparseness (e.g., by anatomical information)\np0::Float64: Prior belief about network sparseness\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.SparseRdcm-NTuple{11, Any}","page":"API reference","title":"RegressionDynamicCausalModeling.SparseRdcm","text":"SparseRdcm(\n    a,\n    c,\n    scans,\n    nr,\n    U,\n    Y,\n    Ep,\n    Conf,\n    hrf,\n    inform_p0,\n    p0\n)\n\n\nConstructor with sanity checks for SparseRdcm. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.SparseRdcm-Tuple{LinearDCM}","page":"API reference","title":"RegressionDynamicCausalModeling.SparseRdcm","text":"SparseRdcm(dcm::LinearDCM; inform_p0, p0) -> SparseRdcm\n\n\nConstruct a SparseRdcm based on a linear DCM.\n\nArguments\n\ndcm: DCM model\ninform_p0::Bool: Inform region specific sparseness (e.g., by anatomical information)\np0::Float64: Prior belief about network sparseness\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.SparseOutput","page":"API reference","title":"RegressionDynamicCausalModeling.SparseOutput","text":"struct SparseOutput <: RegressionDynamicCausalModeling.ModelOutput\n\nOuput after inversion of SparseRdcm.\n\nFields\n\nF::Float64: Negative free energy of the whole model\nF_r::Vector{Float64}: Region specific negative free energy\niter_all::Vector{Int64}: Number of iterations per region until convergence\na_all::Vector{Float64}: Posterior shape parameter for noise\nb_all::Vector{Float64}: Posterior rate parameter for noise\nm_all::Matrix{Float64}: Posterior mean for connectivity parameters\nΣ_all::Vector{SparseArrays.SparseMatrixCSC{Float64, Int64}}: Posterior covariance for connectivity parameters\nz_all::Matrix{Float64}: Posterior for binary indicator variables\ninversion::String: Inversion method\n\n\n\n\n\n","category":"type"},{"location":"api/#Load-and-export-DCMs","page":"API reference","title":"Load and export DCMs","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"load_DCM\nsave_DCM\nexport_to_SPM\nload_example_DCM","category":"page"},{"location":"api/#RegressionDynamicCausalModeling.load_DCM","page":"API reference","title":"RegressionDynamicCausalModeling.load_DCM","text":"load_DCM(path::String;verbose=true)\n\nLoad a DCM from a specified path and return either a linear, bilinear or nonlinear DCM based on the contents of the file.\n\nArguments\n\npath::String: Path to the file to load. Can end with .mat for Matlab files or .jls for data serialized by Julia.\nverbose::Bool: Verbosity\n\nOutput\n\ndcm: DCM struct (can be a linear, bilinear on nonlinear DCM)\n\nExamples\n\njulia> dcm = load_DCM(\"myDCM.jls\";verbose=false)\n\n\n\n\n\n","category":"function"},{"location":"api/#RegressionDynamicCausalModeling.save_DCM","page":"API reference","title":"RegressionDynamicCausalModeling.save_DCM","text":"save_DCM(path, dcm)\n\n\nSave a DCM struct either as .mat or .jls file.\n\nArguments\n\npath::String: Path defining where to save file. If ends with .mat -> saves as Matlab file.   If ends with .jls -> Serializes using Julias Serialization.\ndcm<:DCM: DCM struct (can be a linear, bilinear on nonlinear DCM)\n\nExamples\n\njulia> save_DCM(\"./myDCM.mat\",dcm)\n\n\n\n\n\n","category":"function"},{"location":"api/#RegressionDynamicCausalModeling.export_to_SPM","page":"API reference","title":"RegressionDynamicCausalModeling.export_to_SPM","text":"export_to_SPM(path, rdcm, output)\n\n\nExport a DCM with output after model inversion as an SPM compatible .mat file.\n\nArguments\n\npath::String: Path defining where to save file. Needs to end with .mat\nrdcm::RigidRdcm: An rDCM model.\noutput::RigidOutput: Output after model inversion.\n\nExamples\n\njulia> save_DCM(\"DCM.mat\",rdcm,output)\n\ninfo: Info\nSee SPM compatibility for limitations of this functionality.\n\n\n\n\n\n","category":"function"},{"location":"api/#RegressionDynamicCausalModeling.load_example_DCM","page":"API reference","title":"RegressionDynamicCausalModeling.load_example_DCM","text":"load_example_DCM()\n\n\nLoads an example DCM (linear DCM with 50 regions). The network architecture is based on the S50 structure introduced in Smith et al. (2011). Network modelling methods for FMRI. NeuroImage. This DCM can be used to generate synthetic data.\n\nExamples\n\njulia> dcm = load_example_DCM()\njulia> y_noise, _, _, _ = generate_BOLD(dcm;SNR=10)\n\n\n\n\n\n","category":"function"},{"location":"api/#Data-generation","page":"API reference","title":"Data generation","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"generate_BOLD\npredict","category":"page"},{"location":"api/#RegressionDynamicCausalModeling.generate_BOLD","page":"API reference","title":"RegressionDynamicCausalModeling.generate_BOLD","text":"generate_BOLD(dcm; SNR, TR=NaN, rng=MersenneTwister())\n\nGenerate synthetic BOLD signal timeseries based on a DCM.\n\nArguments\n\ndcm::T where T <: DCM: DCM structure (can be a linear, bilinear on nonlinear DCM)\nSNR::Real: Signal to noise ratio\nTR::Real: Sampling interval in seconds (can be omitted if dt is specified in dcm.Y)\nrng::MersenneTwister: Random number generator for noise sampling.\n\nOutput\n\ny_noise::Matrix{Float64}: BOLD signal timeseries with noise\ny::Matrix{Float64}: BOLD signal timeseries without noise\nx::Matrix{Float64}: Neuronal signal\nh::Vector{Float64}: Hemodynamic response function\n\nExamples\n\njulia> y_noise, _, _, _ = generate_BOLD(load_example_DCM();SNR=10)\n\n\n\n\n\n","category":"function"},{"location":"api/#RegressionDynamicCausalModeling.predict","page":"API reference","title":"RegressionDynamicCausalModeling.predict","text":"predict(\n    rdcm::RDCM,\n    output::RegressionDynamicCausalModeling.ModelOutput\n) -> Any\n\n\nPredict BOLD signal based on parameter estimates after model inversion.\n\nArguments\n\nrdcm: rDCM struct that was used to estimate parameters.\noutput: Output of model inversion.\n\nOutput\n\ny_pred::Matrix{Float64}: predicted BOLD signal\n\n\n\n\n\n","category":"function"},{"location":"api/#Model-inversion","page":"API reference","title":"Model inversion","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"RigidInversionParams\nRigidInversionParams(;::Int64,::Float64)\nSparseInversionParams\nSparseInversionParams(;::Int64,::Float64,::Int64,::Bool)\nOptions\nOptions(::rDCM.AbstractInvParam;::Bool,::Int64,::Bool,::MersenneTwister)\ninvert","category":"page"},{"location":"api/#RegressionDynamicCausalModeling.RigidInversionParams","page":"API reference","title":"RegressionDynamicCausalModeling.RigidInversionParams","text":"struct RigidInversionParams <: RegressionDynamicCausalModeling.AbstractInvParam\n\nSettings for model inversion specific to rigid rDCM (fixed network architecture).\n\nFields\n\nmaxIter::Int64: Maximum number of iterations per region\ntol::Float64: Tolerance for convergence\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.RigidInversionParams-Tuple{}","page":"API reference","title":"RegressionDynamicCausalModeling.RigidInversionParams","text":"RigidInversionParams(; maxIter=500, tol=1.0e-5)\n\nConstructor for RigidInversionParams. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.SparseInversionParams","page":"API reference","title":"RegressionDynamicCausalModeling.SparseInversionParams","text":"struct SparseInversionParams <: RegressionDynamicCausalModeling.AbstractInvParam\n\nSettings for model inversion specific to sparse rDCM.\n\nFields\n\nmaxIter::Int64: Maximum number of iterations per region\ntol::Float64: Tolerance for convergence\nreruns::Int64: Number of reruns\nrestrictInputs::Bool: Whether or not to estimate sparsity for C matrix. If true, the Bernoulli posteriors for the C matrix are not estimated.\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.SparseInversionParams-Tuple{}","page":"API reference","title":"RegressionDynamicCausalModeling.SparseInversionParams","text":"SparseInversionParams(; maxIter=500, tol=1.0e-5, reruns=100, restrictInputs=true)\n\nConstructor for SparseInversionParams. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.Options","page":"API reference","title":"RegressionDynamicCausalModeling.Options","text":"struct Options{T<:RegressionDynamicCausalModeling.AbstractInvParam}\n\nSettings for model inversion.\n\nFields\n\ninvParams::RegressionDynamicCausalModeling.AbstractInvParam: Model specific inversion settings.\nverbose::Int64: Verbosity during inversion (0,1 or 2)\nsynthetic::Bool: Whether or not synthetic data is used.\ntesting::Bool: Used for testing\nrng::Random.MersenneTwister: Random number generator\n\n\n\n\n\n","category":"type"},{"location":"api/#RegressionDynamicCausalModeling.Options-Tuple{RegressionDynamicCausalModeling.AbstractInvParam}","page":"API reference","title":"RegressionDynamicCausalModeling.Options","text":"Options(invParams::T;synthetic::Bool,verbose=1,testing=false,rng=MersenneTwister())\n\nConstructor for Options. See type description for information about arguments.\n\n\n\n\n\n","category":"method"},{"location":"api/#RegressionDynamicCausalModeling.invert","page":"API reference","title":"RegressionDynamicCausalModeling.invert","text":"invert(rdcm::RigidRdcm, opt::Options) -> RigidOutput\n\n\nInvert an rDCM model specified by rdcm with fixed network architecture using specific inversion settings defined in opt.\n\n\n\n\n\ninvert(rdcm::SparseRdcm, opt::Options) -> SparseOutput\n\n\nInvert a sparse rDCM model specified by rdcm using specific inversion settings defined in opt.\n\n\n\n\n\n","category":"function"},{"location":"api/#Index","page":"API reference","title":"Index","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"EditURL = \"../../../examples/Invert_rDCM.jl\"","category":"page"},{"location":"examples/Invert_rDCM/#Inferring-effective-connectivity-with-rDCM","page":"Parameter estimation","title":"Inferring effective connectivity with rDCM","text":"","category":"section"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"In this example we show how to invert a linear DCM with 50 regions using both rDCM with a fixed network architecture and sparse rDCM. In addition to the RegressionDynamicCausalModeling.jl package, Plots.jl and Random.jl are required.","category":"page"},{"location":"examples/Invert_rDCM/#rDCM-with-fixed-network-architecture","page":"Parameter estimation","title":"rDCM with fixed network architecture","text":"","category":"section"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"using RegressionDynamicCausalModeling\nusing Plots\nusing Random\nENV[\"GKSwstype\"] = \"nul\" # operate in \"headless\" mode (for GitHub action) #hide\n\ndcm = load_example_DCM() # this loads the example DCM","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"We can visualize the network architecture of this DCM:","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"heatmap(\n    dcm.a;\n    yflip=true,\n    title=\"Network architecture (A matrix)\",\n    titlefontsize=8,\n    legend=:none,\n    xlabel=\"region from\",\n    ylabel=\"region to\",\n    aspect_ratio=1,\n    size=(350, 350),\n)","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"We can use the DCM to generate synthetic data","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"y_noise, _, _, _ = generate_BOLD(dcm; SNR=3, rng=MersenneTwister(42));\ndcm.Y.y = y_noise # save BOLD signal in DCM\np1 = plot(dcm.Y.y[:, 1]; label=\"true signal\", title=\"BOLD signal of region 1\")","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"Specify options/hyperparameters for inversion","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"opt = Options(RigidInversionParams(); synthetic=true, rng=MersenneTwister(42))\n\nrdcm = RigidRdcm(dcm) # convert the DCM to a rDCM model\nnothing # hide","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"Estimate connectivity","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"output = invert(rdcm, opt)","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"Based on the estimated parameters we can predict the BOLD signal and compare it to the ground truth","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"y_pred = predict(rdcm, output)\nplot(p1, y_pred[:, 1]; label=\"predicted signal\")","category":"page"},{"location":"examples/Invert_rDCM/#sparse-rDCM","page":"Parameter estimation","title":"sparse rDCM","text":"","category":"section"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"We load again the example DCM and generate BOLD signal with an SNR of 3.","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"dcm = load_example_DCM()\ny_noise, _, _, _ = generate_BOLD(dcm; SNR=3, rng=MersenneTwister(42));\ndcm.Y.y = y_noise;\nnothing #hide","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"Convert the linear DCM to a sparse rDCM model with a Bernoulli prior of 0.15 (a priori belief about network sparsity). This means we assume that the a priori probability of a connection being present is low.","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"rdcm = SparseRdcm(dcm; p0=0.15)\nnothing # hide","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"Specify options/hyperparameters for inversion","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"opt = Options(\n    SparseInversionParams(; reruns=100, restrictInputs=true);\n    synthetic=true,\n    rng=MersenneTwister(42),\n)\nnothing # hide","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"Invert the model","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"output = invert(rdcm, opt)","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"We can inspect the Bernoulli posterior for each connection","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"heatmap(\n    output.z_all[:, 1:50];\n    yflip=true,\n    title=\"Posterior over binary indicator variables\",\n    titlefontsize=8,\n    xlabel=\"region from\",\n    ylabel=\"region to\",\n    aspect_ratio=1,\n    size=(350, 350),\n)","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"","category":"page"},{"location":"examples/Invert_rDCM/","page":"Parameter estimation","title":"Parameter estimation","text":"This page was generated using Literate.jl.","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"EditURL = \"../../../examples/Specify_DCM.jl\"","category":"page"},{"location":"examples/Specify_DCM/#How-to-specify-a-DCM","page":"Construct DCM","title":"How to specify a DCM","text":"","category":"section"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"In this toy example we show how to create a simple linear DCM with 3 regions. First we define the network architecture (presence and absence of connections).","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"using RegressionDynamicCausalModeling\nusing Random\n\na = BitMatrix([1 0 1; 1 1 0; 1 0 1])\nc = BitMatrix([0 0; 1 0; 0 1])","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"Number of datapoints per region","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"scans = 300","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"Specify the information about driving inputs","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"off = zeros(100)\non = ones(100)\nu = [off; on; off; on; off]\nu = [u u]\n\nu_dt = 1 / 32 # sampling interval\nname = [\"input1\", \"input2\"]\n\nU = InputU(u, u_dt, name)","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"Next, we specify information about the BOLD signal","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"y = nothing\ny_dt = 0.5 # sampling interval\nname = [\"region1\", \"region2\", \"region3\"]\n\nY = BoldY(y, y_dt, name)","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"The last specification we need are the actual parameter values of the A and C matrix. In this example we generate values by sampling from the prior.","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"Ep = rDCM.TrueParamLinear(a, c; sample=true, rng=MersenneTwister(42))","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"We can now construct our linear DCM and use it for subsequent analysis","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"dcm = LinearDCM(a, c, scans, size(a, 1), U, Y, Ep)","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"","category":"page"},{"location":"examples/Specify_DCM/","page":"Construct DCM","title":"Construct DCM","text":"This page was generated using Literate.jl.","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"EditURL = \"../../../examples/SPM_compat.jl\"","category":"page"},{"location":"examples/SPM_compat/#SPM-compatibility","page":"SPM compatibility","title":"SPM compatibility","text":"","category":"section"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"This package allows to load DCMs that were created by SPM. By loading a DCM.mat file, information that is relevant for rDCM is extracted, all other fields in the DCM struct are ignored. The extracted information is then converted to rDCM's internal representation of a DCM.","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"This internal representation of a DCM can be exported again as a DCM.mat file for e.g., post-processing of inferred rDCM results. However, not all functionality in SPM will be compatible with the exported DCM.mat file because fields that are not relevant for rDCM will be missing. This is mainly due to the fact that some SPM routines expect a B-matrix (i.e., modulatory influences), whereas rDCM is a linear model and contains only A and C parameters. Nevertheless, these incompatibilities relate mainly to, e.g., plotting functionalities in SPM. Based on rDCM results, users will be able to perform the following (and more) analyses:","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"fixed-effects and random-effects Bayesian model selection (BMS)\nfixed-effects and random-effects Bayesian model averaging (BMA)\nBayesian parameter averaging (BPA)\nPost-hoc optimization (search over all reduced models)\nParametric Empirical Bayes (PEB)","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"A DCM.mat file created by SPM can be loaded by","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"dcm = load_DCM(\"DCM.mat\")","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"After model inversion the results can be exported in an SPM compatible format:","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"export_to_SPM(\"DCM.mat\",rdcm,output)","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"","category":"page"},{"location":"examples/SPM_compat/","page":"SPM compatibility","title":"SPM compatibility","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"EditURL = \"https://github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl/blob/main/README.md\"","category":"page"},{"location":"#Regression-Dynamic-Causal-Modeling-(rDCM)","page":"Home","title":"Regression Dynamic Causal Modeling (rDCM)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Dev) (Image: test) (Image: codecov) (Image: Code Style: Blue) (Image: Aqua QA) (Image: JET)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This Julia package implements a variant of dynamic causal modeling (DCM) for fMRI that enables computationally efficient inference on effective (i.e., directed) connectivity parameters among brain regions. Due to its computational efficiency, inversion of large (whole-brain) networks becomes feasible. \nThis package is part of TAPAS which is a collection of software tools developed by the Translational Neuromodeling Unit (TNU) and collaborators. The goal of these tools is to support clinical neuromodeling, particularly computational psychiatry, computational neurology, and computational psychosomatics.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package can be installed using Julia's package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add RegressionDynamicCausalModeling","category":"page"},{"location":"#Minimal-example","page":"Home","title":"Minimal example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following example shows how to simulate synthetic data with an example DCM and then estimate the posterior of the parameters with rDCM. Based on the estimated parameters one can predict the BOLD signal.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using RegressionDynamicCausalModeling # load the package\n\ndcm = load_example_DCM() # load the 50 region example DCM\n\n# Generate synthetic BOLD data with an SNR of 10\ny, _, _, _ = generate_BOLD(dcm; SNR=10)\n\n# Set options for inversion routine\nopt = Options(RigidInversionParams();synthetic=true,verbose=1)\n\n# Convert the linear DCM to a rDCM model with fixed network architecture\nrdcm = RigidRdcm(dcm)\n\n# Invert the model (estimate posterior of parameters)\noutput = invert(rdcm, opt)\n\n# Simulate BOLD signal based on estimated parameters\ny_pred = predict(rdcm,output)","category":"page"},{"location":"#Some-background","page":"Home","title":"Some background","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The regression dynamic causal modeling (rDCM) package implements a variant of DCM for fMRI (Friston et al., 2003) that enables computationally efficient inference on effective (i.e., directed) connectivity among brain regions. This allows rDCM to scale to larger networks and enables whole-brain effective connectivity analyses. \nrDCM was first introduced in Frässle et al. (2017) and then further extended by incorporating sparsity constraints in Frässle et al. (2018). An extension to resting-state fMRI data has been introduced in Frässle et al. (2021). \nImportant note \nThe rDCM framework is in an early stage of development and the method is still subject to limitations. Due to these limitations, the requirements of rDCM in terms of fMRI data quality (i.e., fast repetition time (TR), high signal-to-noise ratio (SNR)) are high - as shown in simulation studies (Frässle et al., 2017; 2018). For data that does not meet these conditions, the method might not give reliable results. It remains the responsibility of the user to ensure that his/her dataset fulfills the requirements.","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This implementation of regression dynamic causal modeling was largely inspired by the original Matlab version.","category":"page"},{"location":"#Citations","page":"Home","title":"Citations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Whenever you use a toolbox from TAPAS in your work, please cite the following paper (main TAPAS reference)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Frässle, S., Aponte, E.A., Bollmann, S., Brodersen, K.H., Do, C.T., Harrison, O.K., Harrison, S.J., Heinzle, J., Iglesias, S., Kasper, L., Lomakina, E.I., Mathys, C., Müller-Schrader, M., Pereira, I., Petzschner, F.H., Raman, S., Schöbi, D., Toussaint, B., Weber, L.A., Yao, Y., Stephan, K.E., 2021. TAPAS: an open-source software package for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry 12, 857.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In addition, please cite the following references if you use the rDCM package:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Frässle, S., Lomakina, E.I., Razi, A., Friston, K.J., Buhmann, J.M., Stephan, K.E., 2017. Regression DCM for fMRI. NeuroImage 155, 406-421.\nFrässle, S., Lomakina, E.I., Kasper, L., Manjaly Z.M., Leff, A., Pruessmann, K.P., Buhmann, J.M., Stephan, K.E., 2018. A generative model of whole-brain effective connectivity. NeuroImage 179, 505-529.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, when using rDCM for resting-state fMRI data, please also cite:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Frässle, S., Harrison, S.J., Heinzle, J., Clementz, B.A., Tamminga, C.A., Sweeney, J.A., Gershon, E.S., Keshavan, M.S., Pearlson, G.D., Powers, A., Stephan, K.E., 2021. Regression dynamic causal modeling for resting-state fMRI. Human Brain Mapping 42, 2159-2180.","category":"page"}]
}
