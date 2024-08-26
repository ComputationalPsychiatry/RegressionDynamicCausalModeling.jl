function test_RigidInversionParams()
    params = RigidInversionParams()

    @test params.maxIter == 500
    @test params.tol == 1.0e-5

    params = RigidInversionParams(;maxIter=100, tol=1.0e-3)

    @test params.maxIter == 100
    @test params.tol == 1.0e-3

    @test_throws ErrorException RigidInversionParams(;maxIter=-3)
    @test_throws ErrorException RigidInversionParams(;tol=-2.0e-3)
end

function test_SparseInversionParams()
    params = SparseInversionParams()

    @test params.maxIter == 500
    @test params.tol == 1.0e-5
    @test params.reruns == 100
    @test params.restrictInputs == true

    @test_throws ErrorException SparseInversionParams(;maxIter=-200)
    @test_throws ErrorException SparseInversionParams(;tol=-1.0e-2)
    @test_throws ErrorException SparseInversionParams(;reruns=-3)
end

function test_InputU()
    U = rDCM.InputU(zeros(100,3),0.5)

    @test all(U.name .== ["u_1","u_2","u_3"])

    U = rDCM.InputU(zeros(100,3),0.5,["allStimuli","a","b"])
    @test U.name[2] == "a"

    @test_throws ErrorException rDCM.InputU(zeros(100,3),-2.5)
    @test_throws ErrorException rDCM.InputU(zeros(100,4),0.5,["a","b","c"])
end

function test_BoldY()
    Y = rDCM.BoldY(zeros(100,3),0.5)

    @test all(Y.name .== ["y_1","y_2","y_3"])

    Y = rDCM.BoldY(zeros(100,3),0.5;name=["reg1","reg2","reg3"])
    @test Y.name[2] == "reg2"

    # test setter function
    @test_throws ErrorException Y.name = ["reg1", "reg2"]
    @test_throws ErrorException Y.name = nothing
    @test_throws ErrorException Y.y = zeros(100,4)
    Y.y = nothing
    @test_throws ErrorException Y.name = ["reg1", "reg2"]
    Y.y = zeros(100,3)
    Y.name = ["a","b","c"]
    @test_throws ErrorException Y.dt = -.5

    @test_throws ErrorException rDCM.BoldY(zeros(100,3),-0.5)
    @test_throws ErrorException rDCM.BoldY(zeros(100,3),0.5;name=["a","b"])

    Y = rDCM.BoldY(zeros(100,3),0.5;name=["y1","y2","y3"])
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
    Tp = rDCM.TrueParamLinear(a,c;rng=MersenneTwister(rDCM.FIXEDSEED))
    A_ref = [-0.5347516797591492 0.0; -1.777533428438784 -0.4983027913744254]
    @test all(Tp.A .== A_ref)

    @test_throws ErrorException rDCM.TrueParamLinear(A,ones(49,3))
    @test_throws ErrorException rDCM.TrueParamLinear(A,C,zeros(50),zeros(49),0.0)
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
    Tp = rDCM.TrueParamBiLinear(a,b,c;rng=MersenneTwister(rDCM.FIXEDSEED))

    A_ref = [-0.6390067190365966 0.0; -0.888766714219392 -0.49321116549770155]
    B_ref = -0.29948409035891055
    C_ref = 1.7778610980573246

    @test all(Tp.A .== A_ref)
    @test Tp.B[1,1,1] == B_ref
    @test Tp.C[1,1] == C_ref

    @test_throws ErrorException rDCM.TrueParamBiLinear(A,B,ones(49,3))
    @test_throws ErrorException rDCM.TrueParamBiLinear(A,B,C,zeros(50),zeros(49),0.0)
end

function test_TrueParamNonLinear()
    A = ones(50,50)
    B = ones(50,50,3)
    C = ones(50,3)
    D = ones(50,50,50)

    @test_throws ErrorException rDCM.TrueParamNonLinear(A,B,ones(49,3),D,zeros(50),zeros(50),0.0)
    @test_throws ErrorException rDCM.TrueParamNonLinear(A,B,C,D,zeros(50),zeros(49),0.0)
end

function test_LinearDCM()
    a = BitArray(ones(2,2))
    a[1,2] = false
    c = BitArray(zeros(2,3))
    c[1,1] = true
    scans = 123
    nr = 2
    nu = size(c,2)
    U = rDCM.InputU(zeros(scans*16,nu),0.03125)
    Y = rDCM.BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamLinear(a,c)

    dcm = LinearDCM(a,c,scans,nr,U,Y,Ep,nothing)

    @test_throws ErrorException dcm.a = BitMatrix(zeros(3,3))
    @test_throws ErrorException dcm.a = BitMatrix(zeros(2,3))
    @test_throws ErrorException dcm.c = BitMatrix(zeros(2,4))
    @test_throws ErrorException dcm.c = BitMatrix(zeros(3,3))
    @test_throws ErrorException dcm.scans = 120
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16,4),0.03125)
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException dcm.Y = rDCM.BoldY(zeros(scans,3),0.5)
    @test_throws ErrorException dcm.Y = rDCM.BoldY(zeros(scans+1,2),0.5)
    @test_throws ErrorException dcm.Ep = rDCM.TrueParamLinear(BitArray(ones(3,3)),BitArray(ones(3,4)))
    @test_throws ErrorException dcm.Ep = rDCM.TrueParamLinear(BitArray(ones(2,2)),BitArray(ones(2,4)))

    @test_throws ErrorException LinearDCM(BitMatrix(ones(3,3)),c,scans,nr,U,Y,Ep,nothing)
    @test_throws ErrorException LinearDCM(a,c,0,nr,U,Y,Ep,nothing)
    conf = rDCM.Confound([1.0, 1.0, 1.0, 1.0],["Constant"])
    @test_throws ErrorException LinearDCM(a,c,scans,nr,U,Y,Ep,conf)

    U_long = rDCM.InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException LinearDCM(a,c,scans,nr,U_long,Y,Ep,nothing)

    # test setter function
    dcm.a = BitMatrix(zeros(2,2))
    dcm.c = BitMatrix(zeros(2,3))
    dcm.scans = 123
    dcm.U = rDCM.InputU(zeros(scans*16,nu),0.03125)
    dcm.Y = rDCM.BoldY(zeros(scans,nr),0.5)
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

    @test_throws ErrorException dcm.a = BitMatrix(zeros(3,3))
    @test_throws ErrorException dcm.a = BitMatrix(zeros(2,3))
    @test_throws ErrorException dcm.b = BitArray(zeros(2,2,4))
    @test_throws ErrorException dcm.b = BitArray(zeros(3,2,3))
    @test_throws ErrorException dcm.c = BitMatrix(zeros(2,4))
    @test_throws ErrorException dcm.c = BitMatrix(zeros(3,3))
    @test_throws ErrorException dcm.scans = 120
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16,4),0.03125)
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException dcm.Y = rDCM.BoldY(zeros(scans,3),0.5)
    @test_throws ErrorException dcm.Y = rDCM.BoldY(zeros(scans+1,2),0.5)
    @test_throws ErrorException dcm.Ep = rDCM.TrueParamBiLinear(BitMatrix(ones(3,3)),BitArray(zeros(3,3,4)),BitMatrix(ones(3,4)))
    @test_throws ErrorException dcm.Ep = rDCM.TrueParamBiLinear(BitMatrix(ones(2,2)),BitArray(zeros(2,2,4)),BitMatrix(ones(2,4)))

    @test_throws ErrorException BiLinearDCM(BitMatrix(ones(3,3)),b,c,scans,nr,U,Y,Ep,nothing)
    @test_throws ErrorException BiLinearDCM(a,b,c,0,nr,U,Y,Ep,nothing)
    @test_throws ErrorException BiLinearDCM(a,b,BitMatrix(zeros(2,4)),scans,nr,U,Y,Ep,nothing)
    conf = rDCM.Confound([1.0, 1.0, 1.0, 1.0],["Constant"])
    @test_throws ErrorException BiLinearDCM(a,b,c,scans,nr,U,Y,Ep,conf)

    U_long = rDCM.InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException BiLinearDCM(a,b,c,scans,nr,U_long,Y,Ep,nothing)

    # test setter function
    dcm.a = BitMatrix(zeros(2,2))
    dcm.b = BitArray(zeros(2,2,3))
    dcm.c = BitMatrix(zeros(2,3))
    dcm.scans = 123
    dcm.U = rDCM.InputU(zeros(scans*16,nu),0.03125)
    dcm.Y = rDCM.BoldY(zeros(scans,nr),0.5)
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
    U = rDCM.InputU(zeros(scans*16,nu),0.03125)
    Y = rDCM.BoldY(zeros(scans,nr),0.5)
    Ep = rDCM.TrueParamNonLinear(a,b,c,d)

    dcm = NonLinearDCM(a,b,c,d,scans,nr,U,Y,Ep,nothing)

    @test_throws ErrorException dcm.a = BitMatrix(zeros(3,3))
    @test_throws ErrorException dcm.a = BitMatrix(zeros(2,3))
    @test_throws ErrorException dcm.b = BitArray(zeros(2,2,4))
    @test_throws ErrorException dcm.b = BitArray(zeros(3,2,3))
    @test_throws ErrorException dcm.c = BitMatrix(zeros(2,4))
    @test_throws ErrorException dcm.c = BitMatrix(zeros(3,3))
    @test_throws ErrorException dcm.d = BitArray(zeros(3,3,3))
    @test_throws ErrorException dcm.d = BitArray(zeros(2,2,3))
    @test_throws ErrorException dcm.scans = 120
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16,4),0.03125)
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException dcm.U = rDCM.InputU(zeros(scans*16,nu),0.03)
    @test_throws ErrorException dcm.Y = rDCM.BoldY(zeros(scans,3),0.5)
    @test_throws ErrorException dcm.Y = rDCM.BoldY(zeros(scans+1,2),0.5)
    @test_throws ErrorException dcm.Ep = rDCM.TrueParamNonLinear(BitMatrix(ones(3,3)),BitArray(zeros(3,3,4)),BitMatrix(ones(3,4)),BitArray(zeros(3,3,3)))
    @test_throws ErrorException dcm.Ep = rDCM.TrueParamNonLinear(BitMatrix(ones(2,2)),BitArray(zeros(2,2,4)),BitMatrix(ones(2,4)),BitArray(zeros(2,2,2)))

    @test_throws ErrorException NonLinearDCM(BitMatrix(ones(3,3)),b,c,d,scans,nr,U,Y,Ep,nothing)
    @test_throws ErrorException NonLinearDCM(a,b,c,d,0,nr,U,Y,Ep,nothing)
    @test_throws ErrorException NonLinearDCM(a,b,BitMatrix(zeros(2,4)),d,scans,nr,U,Y,Ep,nothing)
    conf = rDCM.Confound([1.0, 1.0, 1.0, 1.0],["Constant"])
    @test_throws ErrorException NonLinearDCM(a,b,c,d,scans,nr,U,Y,Ep,conf)

    U_long = rDCM.InputU(zeros(scans*16+1,nu),0.03125)
    @test_throws ErrorException NonLinearDCM(a,b,c,d,scans,nr,U_long,Y,Ep,nothing)

    # test setter function
    dcm.a = BitMatrix(zeros(2,2))
    dcm.b = BitArray(zeros(2,2,3))
    dcm.c = BitMatrix(zeros(2,3))
    dcm.d = BitArray(zeros(2,2,2))
    dcm.scans = 123
    dcm.U = rDCM.InputU(zeros(scans*16,nu),0.03125)
    dcm.Y = rDCM.BoldY(zeros(scans,nr),0.5)
    dcm.Ep = rDCM.TrueParamNonLinear(a,b,c,d)
end

function test_Options()
    @test_throws ErrorException Options(RigidInversionParams();synthetic=true,verbose=-1)
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
    U = rDCM.InputU(zeros(100,2),0.25)
    Ep = rDCM.TrueParamLinear(A,C)
    dcm = LinearDCM(A,C,150,2,U,nothing,Ep)
end

function test_Confound()
    conf = rDCM.Confound([1.0, 1.0, 1.0, 1.0],["Constant"])
    @test all(conf.X0 .== [1.0, 1.0, 1.0, 1.0])
    conf = rDCM.Confound(zeros(10,2))
    @test all(conf.name .== ["Conf_1", "Conf_2"])
end

function test_ModelOutput()
    Σ = [ spzeros(Float64,(2,2)) for _ in 1:4]

    F = 0.0
    F_r = zeros(2)
    iter_all = ones(Int64,2)
    a_all = ones(2)
    b_all = ones(2)
    m_all = zeros(2,2)
    z_all = zeros(2,2) .+ 0.5

    # output of rigid rDCM
    @test_throws ErrorException rDCM.RigidOutput(1.0,F_r,iter_all,      a_all,b_all,  m_all,Σ,"test")
    @test_throws ErrorException rDCM.RigidOutput(F,  F_r,iter_all,      a_all,b_all,  m_all,Σ,"test")
    @test_throws ErrorException rDCM.RigidOutput(F,  F_r,zeros(Int64,2),a_all,b_all,  m_all,Σ,"test")
    @test_throws ErrorException rDCM.RigidOutput(F,  F_r,iter_all,      ones(3),b_all,m_all,Σ,"test")

    # output of sparse rDCM
    @test_throws ErrorException rDCM.SparseOutput(1.0,F_r,iter_all,      a_all,b_all,   m_all,Σ,z_all,          "test")
    @test_throws ErrorException rDCM.SparseOutput(F,  F_r,iter_all,      a_all,zeros(2),m_all,Σ,z_all,          "test")
    @test_throws ErrorException rDCM.SparseOutput(F,  F_r,zeros(Int64,2),a_all,b_all,   m_all,Σ,z_all,          "test")
    @test_throws ErrorException rDCM.SparseOutput(F,  F_r,iter_all,      a_all,b_all,   m_all,Σ,zeros(2,2) .- 1,"test")
    @test_throws ErrorException rDCM.SparseOutput(F,  F_r,iter_all,      a_all,b_all,   m_all,Σ,z_all,          "test")
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
    end
end

test_constructors()
