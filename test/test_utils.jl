
function test_load_DCM(dcm)

    @testset "save/load DCM" begin

        U_name = vec(["u_1","u_2","u_3","u_4","u_5"])
        U_dt = 0.03125
        U_u = zeros((10,7))
        U_u[3:end,1] .= 1.0
        U_u[3:end,end] .= 1.0

        Y_name = vec(["y_1","y_2","y_3","y_4","y_5"])
        Y_dt = 0.5
        Y_y = [-0.0495902545634915	-0.0160198385930071;
        -0.0429054182953994	-0.0135641447362359;
        -0.0370705932010462	-0.0114649325460657;
        -0.0319854078525914	-0.00967398880751166;
        -0.0275602832310189	-0.00814926388392830;
        -0.0237154401097601	-0.00685414022346155;
        -0.0203799883352225	-0.00575702695652323;
        -0.0174911530915810	-0.00483141050923331;
        -0.0149937006836198	-0.00405641035066175;
        -0.0128396105383830	-0.00341773471368428]

        nr = 50
        scans = 2714

        a = BitArray([true	true	false	false	true;
        false	true	true	false	false;
        false	false	true	true	false;
        false	false	false	true	true;
        false	false	false	false	true])

        c = BitArray([0	1	0	0	0;
        0	0	0	0	0;
        0	0	1	0	0;
        0	0	0	0	0;
        0	0	0	0	1])

        @test dcm isa LinearDCM

        @test all(a .== dcm.a[1:5,1:5])
        @test all(c .== dcm.c[1:5,1:5])
        @test scans == dcm.scans
        @test nr    == dcm.nr

        @test all(U_u .== dcm.U.u[1:10,end-6:end])
        @test U_dt == dcm.U.dt
        @test all(cmp.(U_name,dcm.U.name[1:5]) .== 0)

        @test all(Y_y .≈ dcm.Y.y[1:10,1:2])
        @test Y_dt == dcm.Y.dt
        @test all(cmp.(Y_name,dcm.Y.name[1:5]) .== 0)

        @test_throws ErrorException load_DCM(joinpath(rDCM.tmpdir,"DCM_LargeScaleSmith_model1"))

        save_DCM(joinpath(rDCM.tmpdir,"test.jls"),dcm)
        save_DCM(joinpath(rDCM.tmpdir,"test.mat"),dcm)
        @test_throws ErrorException save_DCM(joinpath(rDCM.tmpdir,"test"),dcm)

        dcm_bi = BiLinearDCM(copy(dcm))
        dcm_bi.b[1,1,1] = true
        save_DCM(joinpath(rDCM.tmpdir,"dcm_bi.jls"),dcm_bi)
        save_DCM(joinpath(rDCM.tmpdir,"dcm_bi.mat"),dcm_bi)
        @test load_DCM(joinpath(rDCM.tmpdir,"dcm_bi.jls");verbose=false) isa BiLinearDCM
        @test load_DCM(joinpath(rDCM.tmpdir,"dcm_bi.mat");verbose=false) isa BiLinearDCM

        dcm_non = NonLinearDCM(copy(dcm))
        dcm_non.d[1,1,1] = true
        save_DCM(joinpath(rDCM.tmpdir,"dcm_non.jls"),dcm_non)
        save_DCM(joinpath(rDCM.tmpdir,"dcm_non.mat"),dcm_non)
        @test load_DCM(joinpath(rDCM.tmpdir,"dcm_non.jls");verbose=false) isa NonLinearDCM
        @test load_DCM(joinpath(rDCM.tmpdir,"dcm_non.mat");verbose=false) isa NonLinearDCM
    end
end

function test_print_linearDCM(dcm)
    ref_linear = "Linear DCM
a:     50x50 matrix
c:     50x25 matrix
scans: 2714
nr:    50
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   y:  2714x50 matrix
   dt: 0.5s
   names: y_1,...,y_50
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    C: 50 x 25 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
"

    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", dcm)
    s = String(take!(io))

    @test s == ref_linear

ref_linear = "Linear DCM
a:     50x50 matrix
c:     50x25 matrix
scans: 2714
nr:    50
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   empty
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    C: 50 x 25 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
"

    dcm.Y = nothing
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", dcm)
    s = String(take!(io))

    @test s == ref_linear
end

function test_print_bilinearDCM(dcm)
    ref_bilinear = "Bilinear DCM
a:     50x50 matrix
b:     50x50x25 matrix
c:     50x25 matrix
scans: 2714
nr:    50
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   y: empty
   dt: 0.5s
   names: empty
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    B: 50 x 50 x 25 matrix
    C: 50 x 25 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
"

    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", dcm)
    s = String(take!(io))

    @test s == ref_bilinear

ref_bilinear = "Bilinear DCM
a:     50x50 matrix
b:     50x50x25 matrix
c:     50x25 matrix
scans: 2714
nr:    50
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   empty
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    B: 50 x 50 x 25 matrix
    C: 50 x 25 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
"

    dcm.Y = nothing
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", dcm)
    s = String(take!(io))

    @test s == ref_bilinear
end

function test_print_nonlinearDCM(dcm)
ref_nonlinear = "Nonlinear DCM
a:     50x50 matrix
b:     50x50x25 matrix
c:     50x25 matrix
d:     50x50x50 matrix
scans: 2714
nr:    50
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   y: empty
   dt: 0.5s
   names: empty
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    B: 50 x 50 x 25 matrix
    C: 50 x 25 matrix
    D: 50 x 50 x 50 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
"
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", dcm)
    s = String(take!(io))

    @test s == ref_nonlinear

ref_nonlinear = "Nonlinear DCM
a:     50x50 matrix
b:     50x50x25 matrix
c:     50x25 matrix
d:     50x50x50 matrix
scans: 2714
nr:    50
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   empty
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    B: 50 x 50 x 25 matrix
    C: 50 x 25 matrix
    D: 50 x 50 x 50 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
"

    dcm.Y = nothing
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", dcm)
    s = String(take!(io))

    @test s == ref_nonlinear
end

function test_print_RigidRdcm(rdcm)

ref_rigidrdcm = "rigid rDCM
a:     50x50 matrix
c:     50x26 matrix
scans: 2714
nr:    50
HRF:   43424 element vector
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   y:  2714x50 matrix
   dt: 0.5s
   names: y_1,...,y_50
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    C: 50 x 25 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
---------------------------------------------
Confounds
   X0:    43424x1 matrix
   names: Constant,...,Constant
"

    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", rdcm)
    s = String(take!(io))

    @test s == ref_rigidrdcm
end

function test_print_SparseRdcm(rdcm)

ref_sparserdcm = "sparse rDCM
a:     50x50 matrix
c:     50x26 matrix
scans: 2714
nr:    50
HRF:   43424 element vector
---------------------------------------------
U (Input)
   u:  43424x25 matrix
   dt: 0.03125s
   names: u_1,...,u_25
---------------------------------------------
Y (BOLD signal)
   y:  2714x50 matrix
   dt: 0.5s
   names: y_1,...,y_50
---------------------------------------------
Ep (True parameters)
    A: 50 x 50 matrix
    C: 50 x 25 matrix
    transit: -0.06 ... -0.05
    decay:   0.10 ... 0.10
    epsilon: -0.05
---------------------------------------------
Confounds
   X0:    43424x1 matrix
   names: Constant,...,Constant
---------------------------------------------
inform_p0:false
p0:0.5
"

    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (30, 50)), "text/plain", rdcm)
    s = String(take!(io))

    @test s == ref_sparserdcm
end

function test_spm_compat(dcm)

    # create rDCM struct
    rdcm = RigidRdcm(dcm)

    # set options for inversion
    opt = Options(RigidInversionParams();synthetic=true,verbose=0,testing=true)

    output = invert(rdcm,opt)
    dcm_path = joinpath(rDCM.tmpdir,"dcm_SPM_compat.mat")
    export_to_SPM(dcm_path,rdcm,output)

    dcm2 = load_DCM(dcm_path;verbose=false)

    A_ref = [-0.41548922171677133 0.3949738877798173 0.0;
        0.0 -0.5238266884300515 -0.2711682814395596;
        0.0 0.0 -0.46454372168894786]

    @test all(A_ref .== dcm2.Ep.A[1:3,1:3])

    # test also other fields
    file = matopen(dcm_path)
    DCM_mat = read(file, collect(keys(file))[1])
    close(file)

    # check top level fields
    @test haskey(DCM_mat, "options")
    @test haskey(DCM_mat, "M")
    @test haskey(DCM_mat, "F")
    @test haskey(DCM_mat, "b")
    @test haskey(DCM_mat, "d")
    @test haskey(DCM_mat, "v")
    @test haskey(DCM_mat, "n")
    @test haskey(DCM_mat, "Ce")

    # check lower level fields
    @test haskey(DCM_mat["options"], "nonlinear")
    @test haskey(DCM_mat["options"], "two_state")
    @test haskey(DCM_mat["options"], "stochastic")
    @test haskey(DCM_mat["options"], "centre")
    @test haskey(DCM_mat["options"], "induced")

    @test haskey(DCM_mat["M"], "IS")
    @test haskey(DCM_mat["M"], "pE")
    @test haskey(DCM_mat["M"], "pC")
end

function testUtils()
    @testset verbose=true "Utils" begin

        dcm = load_example_DCM()

        test_load_DCM(copy(dcm))
        @testset "pretty print" begin

            test_print_linearDCM(copy(dcm))
            test_print_bilinearDCM(BiLinearDCM(copy(dcm)))
            test_print_nonlinearDCM(NonLinearDCM(copy(dcm)))
            test_print_RigidRdcm(RigidRdcm(copy(dcm)))
            test_print_SparseRdcm(SparseRdcm(copy(dcm)))
        end

        @testset "SPM compatibility" begin
            test_spm_compat(copy(dcm))
        end
    end
end

testUtils()
