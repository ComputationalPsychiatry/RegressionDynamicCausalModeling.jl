function test_RigidInversionParams()
    params = RigidInversionParams()

    @test params.maxIter == 500
    @test params.tol == 1.0e-5

    params = RigidInversionParams(;maxIter=100, tol=1.0e-3)

    @test params.maxIter == 100
    @test params.tol == 1.0e-3

    # test sanity check of constructor
    @test_throws ErrorException("Invalid parameters.") RigidInversionParams(;maxIter=-3)
    @test_throws ErrorException("Invalid parameters.") RigidInversionParams(;tol=-2.0e-3)
end

function test_SparseInversionParams()
    params = SparseInversionParams()

    @test params.maxIter == 500
    @test params.tol == 1.0e-5
    @test params.reruns == 100
    @test params.restrictInputs == true

    # test sanity check of constructor
    @test_throws ErrorException("Invalid inversion parameters.") SparseInversionParams(;maxIter=-200)
    @test_throws ErrorException("Invalid inversion parameters.") SparseInversionParams(;tol=-1.0e-2)
    @test_throws ErrorException("Invalid inversion parameters.") SparseInversionParams(;reruns=-3)
end

function test_InputU()
    U = InputU(zeros(100,3),0.5)

    @test all(U.name .== ["u_1","u_2","u_3"])

    U = InputU(zeros(100,3),0.5,["allStimuli","a","b"])
    @test U.name[2] == "a"

    # test sanity check of constructor
    @test_throws ErrorException("Sampling rate must be positive.") InputU(zeros(100,3),-2.5)
    @test_throws ErrorException("Size of u and name vector don't match.") InputU(zeros(100,4),0.5,["a","b","c"])
end

function test_BoldY()
    Y = BoldY(zeros(100,3),0.5)

    @test all(Y.name .== ["y_1","y_2","y_3"])

    Y = BoldY(zeros(100,3),0.5;name=["reg1","reg2","reg3"])
    @test Y.name[2] == "reg2"

    # test setter function
    @test_throws ErrorException("Size of y and name vector don't match.") Y.name = ["reg1", "reg2"]
    @test_throws ErrorException("Cannot set name to nothing because y is not nothing.") Y.name = nothing
    @test_throws ErrorException("Size of y and name vector don't match.") Y.y = zeros(100,4)
    Y.y = nothing
    @test_throws ErrorException Y.name = ["reg1", "reg2"]
    Y.y = zeros(100,3)
    Y.name = ["a","b","c"]
    @test_throws ErrorException Y.dt = -.5

    # test sanity check of constructor
    @test_throws ErrorException("Sampling rate must be positive.") BoldY(zeros(100,3),-0.5)
    @test_throws ErrorException("Size of y and name vector don't match.") BoldY(zeros(100,3),0.5;name=["a","b"])

    Y = BoldY(zeros(100,3),0.5;name=["y1","y2","y3"])
    @test all(Y.name .== ["y1", "y2", "y3"])

    # test outer constructor where name is of type Vector{Any}
    rDCM.BoldY(zeros(100,3),0.5;name=[1,"2",3])
end

function test_TrueParamLinear()
    A = ones(50,50)
    C = ones(50,3)
    Tp = rDCM.TrueParamLinear(A,C)

    @test all(Tp.transit .== zeros(50))
    @test all(Tp.decay .== zeros(50))
    @test Tp.epsilon == 0.0

    a = BitArray(ones(2,2))
    a[1,2] = false
    c = BitArray(zeros(2,3))
    Tp = rDCM.TrueParamLinear(a,c;rng=Xoshiro(rDCM.FIXEDSEED))
    A_ref = [-0.45072777489973176 0.0; -3.519434383817597 -0.5546120926381016]
    @test all(Tp.A .== A_ref)

    # test sanity check of constructor
    #@test_throws ErrorException("Size of A and C matrix don't match.") rDCM.TrueParamLinear(A,ones(49,3))
    @test_throws ErrorException("Size of A and C matrix don't match.") rDCM.TrueParamLinear(A,ones(49,3),zeros(50),zeros(50),0.0)
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.TrueParamLinear(A,C,zeros(50),zeros(49),0.0)
end

function test_TrueParamBilinear()
    A = ones(50,50)
    B = ones(50,50,3)
    C = ones(50,3)

    Tp = rDCM.TrueParamBiLinear(A,B,C)
    @test all(Tp.transit .== zeros(50))
    @test all(Tp.decay .== zeros(50))
    @test Tp.epsilon == 0.0

    a = BitArray(ones(2,2))
    a[1,2] = false
    b = BitArray(zeros(2,2,3))
    b[1,1,1] = true
    c = BitArray(zeros(2,3))
    c[1,1] = true
    Tp = rDCM.TrueParamBiLinear(a,b,c;rng=Xoshiro(rDCM.FIXEDSEED))

    A_ref = [-0.30291109959892704 0.0; -1.7597171919087986 -0.7184483705524065]
    B_ref = -0.7332549644348927
    C_ref = -0.7021951987576804


    @test all(Tp.A .== A_ref)
    @test Tp.B[1,1,1] == B_ref
    @test Tp.C[1,1] == C_ref

    # test sanity check of constructor
    #@test_throws ErrorException("Size of A, B or C matrix don't match.") rDCM.TrueParamBiLinear(A,B,ones(49,3)) #TODO: make separate test for priormeanBilinear
    @test_throws ErrorException("Size of A, B or C matrix don't match.") rDCM.TrueParamBiLinear(A,B,ones(49,3),zeros(50),zeros(50),0.0)
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.TrueParamBiLinear(A,B,C,zeros(50),zeros(49),0.0)
end

function test_TrueParamNonLinear()
    A = ones(50,50)
    B = ones(50,50,3)
    C = ones(50,3)
    D = ones(50,50,50)

    # test sanity check of constructor
    @test_throws ErrorException("Size of A, B or C matrix don't match.") rDCM.TrueParamNonLinear(A,B,ones(49,3),D,zeros(50),zeros(50),0.0)
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.TrueParamNonLinear(A,B,C,D,zeros(50),zeros(49),0.0)
end

function test_LinearDCM()
    a = BitArray(ones(2,2))
    a[1,2] = false
    c = BitArray(zeros(2,3))
    c[1,1] = true
    scans = 123
    nr = 2
    nu = size(c,2)
    U = InputU(zeros(scans*16,nu),0.03125)
    Y = BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamLinear(a,c)

    dcm = LinearDCM(a,c,scans,nr,U,Y,Ep)

    @test_throws ErrorException("Number of regions does not match.") dcm.a = BitMatrix(zeros(3,3))
    @test_throws ErrorException("Number of regions does not match.") dcm.a = BitMatrix(zeros(2,3))
    @test_throws ErrorException("Number of inputs does not match.") dcm.c = BitMatrix(zeros(2,4))
    @test_throws ErrorException("Number of regions does not match.") dcm.c = BitMatrix(zeros(3,3))
    @test_throws ErrorException("Number of scans does not match.") dcm.scans = 120
    @test_throws ErrorException("Number of inputs does not match.") dcm.U = InputU(zeros(scans*16,4),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") dcm.U = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") dcm.U = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("Number of regions does not match.") dcm.Y = BoldY(zeros(scans,3),0.5)
    @test_throws ErrorException("Number of scans does not match.") dcm.Y = BoldY(zeros(scans+1,2),0.5)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") dcm.Y = BoldY(zeros(scans,2),0.051)
    @test_throws ErrorException("Number of regions does not match.") dcm.Ep = rDCM.TrueParamLinear(BitArray(ones(3,3)),BitArray(ones(3,4)))
    @test_throws ErrorException("Number of inputs does not match.") dcm.Ep = rDCM.TrueParamLinear(BitArray(ones(2,2)),BitArray(ones(2,4)))

    @test_throws ErrorException("Inconsistent number of regions.") LinearDCM(BitMatrix(ones(3,3)),c,scans,nr,U,Y,Ep,nothing)
    @test_throws ErrorException("Invalid number of scans.") LinearDCM(a,c,0,nr,U,Y,Ep,nothing)
    conf = Confound(ones(scans*16+1),["Constant"])
    @test_throws ErrorException("Confound matrix size and input matrix size don't match.") LinearDCM(a,c,scans,nr,U,Y,Ep,conf)

    Y_long = BoldY(zeros(scans+1,nr),0.5)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") LinearDCM(a,c,scans,nr,U,Y_long,Ep,nothing)

    F_r = zeros(50)
    iter_all = ones(50)
    a_all = ones(50)
    b_all = ones(50)
    m_all = zeros(50,75)
    z_all = zeros(50,75) .+ 0.5
    Σ = [ spzeros(Float64,(75,75)) for _ in 1:50]

    rdcm_rigid = RigidRdcm(dcm)
    out_rigid = rDCM.RigidOutput(0.0,F_r,iter_all,a_all,b_all,m_all,Σ,"test")
    rdcm_sparse = SparseRdcm(dcm;p0=0.5)
    out_sparse = rDCM.SparseOutput(0.0,F_r,iter_all,a_all,b_all,m_all,Σ,z_all,"test")
    @test_throws ErrorException("The output is from a sparse rDCM but the rDCM struct is a rigid rDCM.") LinearDCM(rdcm_rigid,out_sparse)
    @test_throws ErrorException("The output is from a rigid rDCM but the rDCM struct is a sparse rDCM.") LinearDCM(rdcm_sparse,out_rigid)

    # test setter function
    dcm.a = BitMatrix(zeros(2,2))
    dcm.c = BitMatrix(zeros(2,3))
    dcm.scans = 123
    dcm.U = InputU(zeros(scans*16,nu),0.03125)
    dcm.Y = BoldY(zeros(scans,nr),0.5)
    dcm.Ep = rDCM.TrueParamLinear(a,c)
end

function test_BiLinearDCM()
    a = BitMatrix(ones(2,2))
    a[1,2] = false
    b = BitArray(zeros(2,2,3))
    b[1,1,1] = true
    c = BitMatrix(zeros(2,3))
    c[1,1] = true
    scans = 123
    nr = 2
    nu = size(c,2)
    U = rDCM.InputU(zeros(scans*16,nu),0.03125)
    Y = rDCM.BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamBiLinear(a,b,c)

    dcm = BiLinearDCM(a,b,c,scans,nr,U,Y,Ep,nothing)

    @test_throws ErrorException("Number of regions does not match.") dcm.a = BitMatrix(zeros(3,3))
    @test_throws ErrorException("Number of regions does not match.") dcm.a = BitMatrix(zeros(2,3))
    @test_throws ErrorException("Number of inputs does not match.") dcm.b = BitArray(zeros(2,2,4))
    @test_throws ErrorException("Number of regions does not match.") dcm.b = BitArray(zeros(3,2,3))
    @test_throws ErrorException("Number of inputs does not match.") dcm.c = BitMatrix(zeros(2,4))
    @test_throws ErrorException("Number of regions does not match.") dcm.c = BitMatrix(zeros(3,3))
    @test_throws ErrorException("Number of scans does not match.") dcm.scans = 120
    @test_throws ErrorException("Number of inputs does not match.") dcm.U = InputU(zeros(scans*16,4),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") dcm.U = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") dcm.U = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("Number of regions does not match.") dcm.Y = BoldY(zeros(scans,3),0.5)
    @test_throws ErrorException("Number of scans does not match.") dcm.Y = BoldY(zeros(scans+1,2),0.5)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") dcm.Y = BoldY(zeros(scans,2),0.51)
    @test_throws ErrorException("Number of regions does not match.") dcm.Ep = rDCM.TrueParamBiLinear(BitMatrix(ones(3,3)),BitArray(zeros(3,3,4)),BitMatrix(ones(3,4)))
    @test_throws ErrorException("Number of inputs does not match.") dcm.Ep = rDCM.TrueParamBiLinear(BitMatrix(ones(2,2)),BitArray(zeros(2,2,4)),BitMatrix(ones(2,4)))

    @test_throws ErrorException("Inconsistent number of regions.") BiLinearDCM(BitMatrix(ones(3,3)),b,c,scans,nr,U,Y,Ep,nothing)
    @test_throws ErrorException("Invalid number of scans.") BiLinearDCM(a,b,c,0,nr,U,Y,Ep,nothing)
    @test_throws ErrorException("Number of inputs don't match.") BiLinearDCM(a,b,BitMatrix(zeros(2,4)),scans,nr,U,Y,Ep,nothing)
    conf = rDCM.Confound(zeros(scans*16+1),["Constant"])
    @test_throws ErrorException("Confound matrix size and input matrix size don't match.") BiLinearDCM(a,b,c,scans,nr,U,Y,Ep,conf)

    U_long = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") BiLinearDCM(a,b,c,scans,nr,U_long,Y,Ep,nothing)

    # wrong sampling rate
    U_wrong_dt = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") BiLinearDCM(a,b,c,scans,nr,U_wrong_dt,Y,Ep,nothing)

    # test setter function
    dcm.a = BitMatrix(zeros(2,2))
    dcm.b = BitArray(zeros(2,2,3))
    dcm.c = BitMatrix(zeros(2,3))
    dcm.scans = 123
    dcm.U = InputU(zeros(scans*16,nu),0.03125)
    dcm.Y = BoldY(zeros(scans,nr),0.5)
    dcm.Ep = rDCM.TrueParamBiLinear(a,b,c)
end

function test_NonLinearDCM()
    a = BitMatrix(ones(2,2))
    a[1,2] = false
    b = BitArray(zeros(2,2,3))
    b[1,1,1] = true
    c = BitMatrix(zeros(2,3))
    c[1,1] = true
    d = BitArray(zeros(2,2,2))
    d[1,1,1] = true
    scans = 123
    nr = 2
    nu = size(c,2)
    U = InputU(zeros(scans*16,nu),0.03125)
    Y = BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamNonLinear(a,b,c,d)

    dcm = NonLinearDCM(a,b,c,d,scans,nr,U,Y,Ep,nothing)

    @test_throws ErrorException("Number of regions does not match.") dcm.a = BitMatrix(zeros(3,3))
    @test_throws ErrorException("Number of regions does not match.") dcm.a = BitMatrix(zeros(2,3))
    @test_throws ErrorException("Number of inputs does not match.") dcm.b = BitArray(zeros(2,2,4))
    @test_throws ErrorException("Number of regions does not match.") dcm.b = BitArray(zeros(3,2,3))
    @test_throws ErrorException("Number of inputs does not match.") dcm.c = BitMatrix(zeros(2,4))
    @test_throws ErrorException("Number of regions does not match.") dcm.c = BitMatrix(zeros(3,3))
    @test_throws ErrorException("Number of regions does not match.") dcm.d = BitArray(zeros(3,3,3))
    @test_throws ErrorException("Number of regions does not match.") dcm.d = BitArray(zeros(2,2,3))
    @test_throws ErrorException("Number of scans does not match.") dcm.scans = 120
    @test_throws ErrorException("Number of inputs does not match.") dcm.U = InputU(zeros(scans*16,4),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") dcm.U = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") dcm.U = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("Number of regions does not match.") dcm.Y = BoldY(zeros(scans,3),0.5)
    @test_throws ErrorException("Number of scans does not match.") dcm.Y = BoldY(zeros(scans+1,2),0.5)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") dcm.Y = BoldY(zeros(scans,2),0.51)
    @test_throws ErrorException("Number of regions does not match.") dcm.Ep = rDCM.TrueParamNonLinear(BitMatrix(ones(3,3)),BitArray(zeros(3,3,4)),BitMatrix(ones(3,4)),BitArray(zeros(3,3,3)))
    @test_throws ErrorException("Number of inputs does not match.") dcm.Ep = rDCM.TrueParamNonLinear(BitMatrix(ones(2,2)),BitArray(zeros(2,2,4)),BitMatrix(ones(2,4)),BitArray(zeros(2,2,2)))

    @test_throws ErrorException("Inconsistent number of regions.") NonLinearDCM(BitMatrix(ones(3,3)),b,c,d,scans,nr,U,Y,Ep,nothing)
    @test_throws ErrorException("Invalid number of scans.") NonLinearDCM(a,b,c,d,0,nr,U,Y,Ep,nothing)
    @test_throws ErrorException("Number of inputs don't match.") NonLinearDCM(a,b,BitMatrix(zeros(2,4)),d,scans,nr,U,Y,Ep,nothing)
    conf = Confound(zeros(scans*16+1),["Constant"])
    @test_throws ErrorException("Confound matrix size and input matrix size don't match.") NonLinearDCM(a,b,c,d,scans,nr,U,Y,Ep,conf)

    U_long = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") NonLinearDCM(a,b,c,d,scans,nr,U_long,Y,Ep,nothing)

    # wrong sampling rate
    U_wrong_dt = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") NonLinearDCM(a,b,c,d,scans,nr,U_wrong_dt,Y,Ep,nothing)

    # test setter function
    dcm.a = BitMatrix(zeros(2,2))
    dcm.b = BitArray(zeros(2,2,3))
    dcm.c = BitMatrix(zeros(2,3))
    dcm.d = BitArray(zeros(2,2,2))
    dcm.scans = 123
    dcm.U = InputU(zeros(scans*16,nu),0.03125)
    dcm.Y = BoldY(zeros(scans,nr),0.5)
    dcm.Ep = rDCM.TrueParamNonLinear(a,b,c,d)
end

function test_Options()
    @test_throws ErrorException("Verbosity level needs to be an integer between 0 and 2.") Options(RigidInversionParams();synthetic=true,verbose=-1)
end

function test_TrueParam()
    a = BitMatrix(zeros(2,2))
    a[1,1] = true
    c = BitMatrix(zeros(2,2))
    c[1,1] = true
    Ep = rDCM.TrueParamLinear(a,c;sample=false)

    @test Ep.A[1,1] == -.5
    @test Ep.C[1,1] == 0

    b = BitArray(zeros(2,2,2))
    b[1,1,1] = true
    Ep = rDCM.TrueParamLinear(rDCM.TrueParamBiLinear(a,b,c;sample=false))
    @test Ep.A[1,1] == -.5
    @test Ep.C[1,1] == 0

    d = BitArray(zeros(2,2,2))
    d[1,1,1] = true
    Ep = rDCM.TrueParamLinear(rDCM.TrueParamNonLinear(a,b,c,d;sample=false))
    @test Ep.A[1,1] == -.5
    @test Ep.C[1,1] == 0

    Ep = rDCM.TrueParamNonLinear(rDCM.TrueParamBiLinear(a,b,c;sample=false))
    @test all(Ep.D .== zeros(2,2,2))

    Ep = rDCM.TrueParamNonLinear(Ep.A,Ep.C)
    @test all(Ep.B .== zeros(2,2,2))
    @test all(Ep.D .== zeros(2,2,2))
end

function test_DCM()
    dcm = load_example_DCM()
    B = zeros(50,50,size(dcm.c,2))
    B[1,1,1] = 0.5
    dcm_bi = copy(BiLinearDCM(dcm,B))
    @test dcm_bi.b[1,1,1] == true

    dcm_bi = rDCM.load_DCM_bilinear()
    dcm_non = rDCM.load_DCM_nonlinear()

    dcm = LinearDCM(dcm_bi)
    dcm = LinearDCM(dcm_non)
    dcm_non = NonLinearDCM(dcm_bi)
    A = zeros(2,2)
    A[1,1] = 0.8
    C = zeros(2,2)
    C[1,1] = -.4
    U = InputU(zeros(100,2),0.25)
    Ep = rDCM.TrueParamLinear(A,C)
    dcm = LinearDCM(A,C,150,2,U,nothing,Ep)
end

function test_Confound()
    conf = Confound([1.0, 1.0, 1.0, 1.0],["Constant"])
    @test all(conf.X0 .== [1.0, 1.0, 1.0, 1.0])
    conf = Confound(zeros(10,2))
    @test all(conf.name .== ["Conf_1", "Conf_2"])
    @test_throws ErrorException("Size of X0 and name vector don't match.") Confound(zeros(3,2),["Const"])
end

function test_ModelOutput()
    Σ = [ spzeros(Float64,(2,2)) for _ in 1:2]

    F = 0.0
    F_r = zeros(2)
    iter_all = ones(Int,2)
    a_all = ones(2)
    b_all = ones(2)
    m_all = ones(2,2)
    z_all = zeros(2,2) .+ 0.5

    # output of rigid rDCM
    @test_throws ErrorException("Sum of region-wise neg. free energies don't sum up to overall neg. free energy.") rDCM.RigidOutput(1.0,F_r,iter_all,      a_all,b_all,   m_all,Σ,"test")
    @test_throws ErrorException("Found invalid values of the posterior Gamma distribution.") rDCM.RigidOutput(F,  F_r,iter_all,      a_all,zeros(2),m_all,Σ,"test")
    @test_throws ErrorException("Invalid number of iterations") rDCM.RigidOutput(F,  F_r,zeros(Int,2),a_all,b_all,   m_all,Σ,"test")
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.RigidOutput(F,  F_r,iter_all,      ones(3),b_all, m_all,Σ,"test")
    @test_throws ErrorException("One or more covariance matrices are not positive definite.") rDCM.RigidOutput(F, F_r, iter_all, a_all,b_all, m_all, [sparse([1 2; 3 4]) for _ in 1:2],"test")

    # output of sparse rDCM
    @test_throws ErrorException("Sum of region-wise neg. free energies don't sum up to overall neg. free energy.") rDCM.SparseOutput(1.0,F_r,iter_all,      a_all,b_all,   m_all,Σ,z_all,          "test")
    @test_throws ErrorException("Found invalid values of the posterior Gamma distribution.") rDCM.SparseOutput(F,  F_r,iter_all,      a_all,zeros(2),m_all,Σ,z_all,          "test")
    @test_throws ErrorException("Invalid number of iterations") rDCM.SparseOutput(F,  F_r,zeros(Int,2),a_all,b_all,   m_all,Σ,z_all,          "test")
    @test_throws ErrorException("Invalid probabilities in posterior Bernoulli.") rDCM.SparseOutput(F,  F_r,iter_all,      a_all,b_all,   m_all,Σ,zeros(2,2) .- 1,"test")
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.SparseOutput(F,  F_r,iter_all,      a_all,b_all,   m_all,[ spzeros(Float64,(2,2)) for _ in 1:3],z_all,          "test")
    @test_throws ErrorException("One or more covariance matrices are not positive definite.") rDCM.SparseOutput(F, F_r, iter_all, a_all, b_all, m_all, [sparse([1 2;3 4]) for _ in 1:2], z_all, "test")
end

function test_RigiRdcm()
    a = BitArray(ones(2,2))
    a[1,2] = false
    c = BitArray(zeros(2,3))
    c[1,1] = true
    scans = 123
    nr = 2
    nu = size(c,2)
    U = InputU(zeros(scans*16,nu),0.03125)
    Y = BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamLinear(a,c)
    conf = Confound(ones(scans*16),["Constant"])
    hrf = zeros(scans*16)

    @test_throws ErrorException("Dimension mismatch.") RigidRdcm(BitMatrix(ones(3,3)),c,scans,nr,U,Y,Ep,conf,hrf)
    @test_throws ErrorException("Number of inputs don't match.") RigidRdcm(a,BitMatrix(ones(2,4)),scans,nr,U,Y,Ep,conf,hrf)

    # wrong sampling rate of U
    U_wrong_dt = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") RigidRdcm(a,c,scans,nr,U_wrong_dt,Y,Ep,conf,hrf)
    # wrong length of U
    U_long = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") RigidRdcm(a,c,scans,nr,U_long,Y,Ep,conf,hrf)

    # test outer constructors
    dcm = LinearDCM(a,c,scans,nr,U,Y,Ep,conf)
    RigidRdcm(dcm)
    # U is nothing
    dcm = LinearDCM(a,c,scans,nr,nothing,Y,Ep,conf)
    RigidRdcm(dcm)
    # U and conf are nothing
    dcm = LinearDCM(a,c,scans,nr,nothing,Y,Ep,nothing)
    RigidRdcm(dcm)

end

function test_SparseRdcm()
    a = BitArray(ones(2,2))
    a[1,2] = false
    c = BitArray(zeros(2,3))
    c[1,1] = true
    scans = 123
    nr = 2
    nu = size(c,2)
    U = InputU(zeros(scans*16,nu),0.03125)
    Y = BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamLinear(a,c)
    conf = Confound(ones(3),["Constant"])
    hrf = zeros(scans*16)
    p0 = 0.5

    @test_throws ErrorException("Dimension mismatch.") SparseRdcm(BitMatrix(ones(3,3)),c,scans,nr,U,Y,Ep,conf,hrf,true,p0)
    @test_throws ErrorException("Number of inputs don't match.") SparseRdcm(a,BitMatrix(ones(2,4)),scans,nr,U,Y,Ep,conf,hrf,true,p0)
    @test_throws ErrorException("p0 is not a proper Bernoulli parameter.") SparseRdcm(a,c,scans,nr,U,Y,Ep,conf,hrf,true,1.5)
    U_wrong_dt = InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate of the input U (u_dt). Cannot proceed.") SparseRdcm(a,c,scans,nr,U_wrong_dt,Y,Ep,conf,hrf,true,p0)
    # wrong length of U
    U_long = InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException("Length of BOLD signal and driving input u is inconsistent.") SparseRdcm(a,c,scans,nr,U_long,Y,Ep,conf,hrf,true,p0)

end

function test_priorMean()
    nr = 3
    nu = 4
    B = zeros(nr,nr,nu)
    C = zeros(nr,nu)
    D = zeros(nr,nr,nr)
    transit = zeros(nr)
    decay = zeros(nr)
    epsilon = 0.0
    A_wrong = zeros(nr+1,nr+1)
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.priorMeanLinear(A_wrong,C,transit,decay,epsilon)
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.priorMeanBiLinear(A_wrong,B,C,transit,decay,epsilon)
    @test_throws ErrorException("Inconsistent number of regions.") rDCM.priorMeanNonLinear(A_wrong,B,C,D,transit,decay,epsilon)
end

function test_constructors()
    @testset verbose=true "Constructors" begin
        test_RigidInversionParams()
        test_SparseInversionParams()
        test_InputU()
        test_BoldY()
        test_TrueParamLinear()
        test_TrueParamBilinear()
        test_TrueParamNonLinear()
        test_LinearDCM()
        test_BiLinearDCM()
        test_NonLinearDCM()
        test_Options()
        test_TrueParam()
        test_DCM()
        test_Confound()
        test_ModelOutput()
        test_RigiRdcm()
        test_SparseRdcm()
        test_priorMean()
    end
end

test_constructors()
