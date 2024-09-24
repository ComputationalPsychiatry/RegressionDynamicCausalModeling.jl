# # How to specify a DCM/rDCM
# ## Simple DCM
# In this toy example we show how to create a simple linear DCM with 3 regions.
# First we define the network architecture (presence and absence of connections).
using RegressionDynamicCausalModeling
using Random

a = BitMatrix([1 0 1; 1 1 0; 1 0 1])
c = BitMatrix([0 0; 1 0; 0 1])

# Number of datapoints per region
scans = 300

# Specify the information about driving inputs
off = zeros(100)
on = ones(100)
u = [off; on; off; on; off]
u = [u u]

u_dt = 1 / 32 # sampling interval
name = ["input1", "input2"]

U = InputU(u, u_dt, name)

# Next, we specify information about the BOLD signal
y = nothing
y_dt = 0.5 # sampling interval
name = ["region1", "region2", "region3"]

Y = BoldY(y, y_dt, name)

# The last specification we need are the actual parameter values of the A and C matrix.
# In this example we generate values by sampling from the prior.
Ep = rDCM.TrueParamLinear(a, c; sample=true, rng=MersenneTwister(42))

# We can now construct our linear DCM and use it for subsequent analysis
dcm = LinearDCM(a, c, scans, size(a, 1), U, Y, Ep)

# ## Confounds
# It is possible to specify confounds. The following example shows how this can be done.
X0 = zeros(500,3) # specify confound matrix (fist dimension is time, second dimension is number of confounds)
conf_names = ["conf1", "conf2", "conf3"]
conf = Confound(X0,conf_names)

# Construct DCM and convert to rDCM structure
dcm = LinearDCM(a, c, scans, size(a, 1), U, Y, Ep,conf)
rdcm = RigidRdcm(dcm)
