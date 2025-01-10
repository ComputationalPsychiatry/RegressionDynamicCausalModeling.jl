# # SPM compatibility
#
# This package allows to load DCMs that were created by [SPM](https://www.fil.ion.ucl.ac.uk/spm/) [^1].
# By loading a `DCM.mat` file, information that is relevant for rDCM is extracted, all other
# fields in the DCM *struct* are ignored. The extracted information is then converted to
# rDCM's internal representation of a DCM.

# This internal representation of a DCM can be exported again as a `DCM.mat` file for e.g.,
# post-processing of inferred rDCM results. However, not all functionality in SPM will be
# compatible with the exported `DCM.mat` file because fields that are not relevant for rDCM
# will be missing. This is mainly due to the fact that some SPM routines expect a B-matrix
# (i.e., modulatory influences), whereas rDCM is a linear model and contains only A and C
# parameters. Nevertheless, these incompatibilities relate mainly to, e.g., plotting
# functionalities in SPM.
# Based on rDCM results, users will be able to perform the following (and more) analyses:
# - fixed-effects and random-effects Bayesian model selection (BMS)
# - fixed-effects and random-effects Bayesian model averaging (BMA)
# - Bayesian parameter averaging (BPA)
# - Post-hoc optimization (search over all reduced models)
# - Parametric Empirical Bayes (PEB)


# A `DCM.mat` file created by SPM can be loaded by
# ```
# dcm = load_DCM("DCM.mat")
# ```

# After model inversion the results can be exported in an SPM compatible format:
# ```
# export_to_SPM("DCM.mat",rdcm,output)
# ```

# [^1]: tested with SPM12 (7771)
