# API reference

```@docs
RegressionDynamicCausalModeling
```

## DCMs
```@docs
DCM
LinearDCM
LinearDCM(a,c,scans,nr,U,Y,Ep,Conf)
LinearDCM(::Matrix{Number},::Matrix{Number},::Int64,::Int64,::InputU,::Union{BoldY,Nothing},::rDCM.TrueParamLinear)
LinearDCM(::DCM)
LinearDCM(::RDCM,::rDCM.ModelOutput)
BiLinearDCM
BiLinearDCM(a,b,c,scans,nr,U,Y,Ep,Conf)
NonLinearDCM
NonLinearDCM(a,b,c,d,scans,nr,U,Y,Ep,Conf)
InputU
BoldY
Confound
```

## rDCMs
```@docs
RDCM
RigidRdcm
RigidRdcm(a,c,scans,nr,U,Y,Ep,Conf,hrf)
RigidRdcm(::LinearDCM)
RigidOutput
SparseRdcm
SparseRdcm(a,c,scans,nr,U,Y,Ep,Conf,hrf,inform_p0,p0)
SparseRdcm(::LinearDCM;::Bool,::Float64)
SparseOutput
```

## Load and export DCMs
```@docs
load_DCM
save_DCM
export_to_SPM
load_example_DCM
```

## Data generation
```@docs
generate_BOLD
predict
```

## Model inversion
```@docs
RigidInversionParams
RigidInversionParams(;::Int64,::Float64)
SparseInversionParams
SparseInversionParams(;::Int64,::Float64,::Int64,::Bool)
Options
Options(::rDCM.AbstractInvParam;::Bool,::Int64,::Bool,::MersenneTwister)
invert
```

## Index
```@index
```
