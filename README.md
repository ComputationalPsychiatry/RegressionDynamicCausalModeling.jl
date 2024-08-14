# Regression Dynamic Causal Modeling (rDCM)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ComputationalPsychiatry.github.io/RegressionDynamicCausalModeling.jl/dev/)
[![test](https://github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl/actions/workflows/test.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl/graph/badge.svg?token=6GHXMZTQEJ)](https://codecov.io/gh/ComputationalPsychiatry/RegressionDynamicCausalModeling.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%E2%9C%88%EF%B8%8F%20tested%20with%20-%20JET.jl%20-%20red)](https://github.com/aviatesk/JET.jl)
[![RegressionDynamicCausalModeling Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FRegressionDynamicCausalModeling&query=total_requests&suffix=%2Fmonth&label=Downloads)](http://juliapkgstats.com/pkg/RegressionDynamicCausalModeling)

This Julia package implements a variant of dynamic causal modeling (DCM) for fMRI that enables computationally efficient inference on effective (i.e., directed) connectivity parameters among brain regions. Due to its computational efficiency, inversion of large (whole-brain) networks becomes feasible. \
This package is part of [TAPAS](https://translationalneuromodeling.github.io/tapas/) which is a collection of software tools developed by the [Translational Neuromodeling Unit](https://www.tnu.ethz.ch/en/home) (TNU) and collaborators. The goal of these tools is to support clinical neuromodeling, particularly computational psychiatry, computational neurology, and computational psychosomatics.

## Getting started

This package can be installed using Julia's package manager:
```julia
pkg> add RegressionDynamicCausalModeling
```

### Minimal example
The following example shows how to simulate synthetic data with an example DCM and then estimate the posterior of the parameters with rDCM. Based on the estimated parameters one can predict the BOLD signal.


```julia
using RegressionDynamicCausalModeling # load the package

dcm = load_example_DCM() # load the 50 region example DCM

# Generate synthetic BOLD data with an SNR of 10
y, _, _, _ = generate_BOLD(dcm; SNR=10)

# Set options for inversion routine
opt = Options(RigidInversionParams();synthetic=true,verbose=1)

# Convert the linear DCM to a rDCM model with fixed network architecture
rdcm = RigidRdcm(dcm)

# Invert the model (estimate posterior of parameters)
output = invert(rdcm, opt)

# Simulate BOLD signal based on estimated parameters
y_pred = predict(rdcm,output)
```

Detailed documentation can be found [here](https://ComputationalPsychiatry.github.io/RegressionDynamicCausalModeling.jl/dev/).

## Some background
The regression dynamic causal modeling (rDCM) package implements a variant of DCM for fMRI ([Friston et al., 2003](https://pubmed.ncbi.nlm.nih.gov/12948688/)) that enables computationally efficient inference on effective (i.e., directed) connectivity among brain regions. This allows rDCM to scale to larger networks and enables whole-brain effective connectivity analyses. \
rDCM was first introduced in [Frässle et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28259780/) and then further extended by incorporating sparsity constraints in [Frässle et al. (2018)](https://www.sciencedirect.com/science/article/pii/S1053811918304762). An extension to resting-state fMRI data has been introduced in [Frässle et al. (2021)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8046067/). \
**Important note** \
The rDCM framework is in an early stage of development and the method is still subject to limitations. Due to these limitations, the requirements of rDCM in terms of fMRI data quality (i.e., fast repetition time (TR), high signal-to-noise ratio (SNR)) are high - as shown in simulation studies (Frässle et al., 2017; 2018). For data that does not meet these conditions, the method might not give reliable results. It remains the responsibility of the user to ensure that his/her dataset fulfills the requirements.

## Acknowledgements
This implementation of regression dynamic causal modeling was largely inspired by the original [Matlab version](https://github.com/translationalneuromodeling/tapas/tree/master/rDCM).

## Citations
Whenever you use a toolbox from [TAPAS](https://translationalneuromodeling.github.io/tapas/) in your work, please cite the following paper (main TAPAS reference)

- Frässle, S., Aponte, E.A., Bollmann, S., Brodersen, K.H., Do, C.T., Harrison, O.K., Harrison, S.J., Heinzle, J., Iglesias, S., Kasper, L., Lomakina, E.I., Mathys, C., Müller-Schrader, M., Pereira, I., Petzschner, F.H., Raman, S., Schöbi, D., Toussaint, B., Weber, L.A., Yao, Y., Stephan, K.E., 2021. TAPAS: an open-source software package for Translational Neuromodeling and Computational Psychiatry. *Frontiers in Psychiatry* 12, 857.

In addition, please cite the following references if you use the rDCM package:
- Frässle, S., Lomakina, E.I., Razi, A., Friston, K.J., Buhmann, J.M., Stephan, K.E., 2017. Regression DCM for fMRI. NeuroImage 155, 406-421.
- Frässle, S., Lomakina, E.I., Kasper, L., Manjaly Z.M., Leff, A., Pruessmann, K.P., Buhmann, J.M., Stephan, K.E., 2018. A generative model of whole-brain effective connectivity. *NeuroImage* 179, 505-529.

Finally, when using rDCM for resting-state fMRI data, please also cite:
- Frässle, S., Harrison, S.J., Heinzle, J., Clementz, B.A., Tamminga, C.A., Sweeney, J.A., Gershon, E.S., Keshavan, M.S., Pearlson, G.D., Powers, A., Stephan, K.E., 2021. Regression dynamic causal modeling for resting-state fMRI. *Human Brain Mapping* 42, 2159-2180.