# # Inferring effective connectivity with rDCM

# In this example we show how to invert a linear DCM with 50 regions using both rDCM with a
# fixed network architecture and sparse rDCM. In addition to the
# [RegressionDynamicCausalModeling.jl](https://github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl)
# package, [Plots.jl](https://docs.juliaplots.org/stable/)
# and [Random.jl](https://docs.julialang.org/en/v1/stdlib/Random/) are required.

# ## rDCM with fixed network architecture
using RegressionDynamicCausalModeling
using Plots
using Random
ENV["GKSwstype"] = "nul" # operate in "headless" mode (for GitHub action) #hide

dcm = load_example_DCM() # this loads the example DCM

# We can visualize the network architecture of this DCM:
heatmap(
    dcm.a;
    yflip=true,
    title="Network architecture (A matrix)",
    titlefontsize=8,
    legend=:none,
    xlabel="region from",
    ylabel="region to",
    aspect_ratio=1,
    size=(350, 350),
)

# We can use the DCM to generate synthetic data
y_noise, _, _, _ = generate_BOLD(dcm; SNR=3, rng=MersenneTwister(42));
dcm.Y.y = y_noise # save BOLD signal in DCM
p1 = plot(dcm.Y.y[:, 1]; label="true signal", title="BOLD signal of region 1")

# Specify options/hyperparameters for inversion
opt = Options(RigidInversionParams(); synthetic=true, rng=MersenneTwister(42))

rdcm = RigidRdcm(dcm) # convert the DCM to a rDCM model
nothing # hide

# Estimate connectivity
output = invert(rdcm, opt)

# Based on the estimated parameters we can predict the BOLD signal and compare it to the
# ground truth
y_pred = predict(rdcm, output)
plot(p1, y_pred[:, 1]; label="predicted signal")

# ## sparse rDCM
# We load again the example DCM and generate BOLD signal with an SNR of 3.
dcm = load_example_DCM()
y_noise, _, _, _ = generate_BOLD(dcm; SNR=3, rng=MersenneTwister(42));
dcm.Y.y = y_noise;

# Convert the linear DCM to a sparse rDCM model with a Bernoulli prior of 0.15
# (a priori belief about network sparsity). This means we assume that the a priori
# probability of a connection being present is low.
rdcm = SparseRdcm(dcm; p0=0.15)
nothing # hide

# Specify options/hyperparameters for inversion
opt = Options(
    SparseInversionParams(; reruns=100, restrictInputs=true);
    synthetic=true,
    rng=MersenneTwister(42),
)
nothing # hide

# Invert the model
output = invert(rdcm, opt)

# We can inspect the Bernoulli posterior for each connection
heatmap(
    output.z_all[:, 1:50];
    yflip=true,
    title="Posterior over binary indicator variables",
    titlefontsize=8,
    xlabel="region from",
    ylabel="region to",
    aspect_ratio=1,
    size=(350, 350),
)
