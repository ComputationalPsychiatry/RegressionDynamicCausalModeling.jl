
function test_HRF(dcm)

    @testset "Create HRF" begin
        N = 43424
        u_dt = 1/32
        hrf = rDCM.get_hrf(N,u_dt)

        hrf_ref = [0,
        0,
        0,
        0,
        -3.29103402287245e-07,
        -1.11865515089661e-06,
        -2.32331983047075e-06,
        -3.74068920321818e-06,
        -5.02950335036502e-06,
        -5.72656552901520e-06,
        -5.26242693096037e-06,
        -2.97591369046443e-06,
        1.87243514211530e-06,
        1.00879536256134e-05,
        2.25303971330018e-05,
        4.01033944357754e-05,
        6.37447482576026e-05,
        9.44175320615305e-05,
        0.000133101933616547,
        0.000180787798599450]

        @test all(hrf_ref .≈ hrf[1:20])

        u_dt = 1/64
        hrf = rDCM.get_hrf(N,u_dt)

        hrf_ref = [0,
        0,
        0,
        0,
        -2.05702753665094e-08,
        -7.61047937254489e-08,
        -1.75151507901461e-07,
        -3.20731039526407e-07,
        -5.10644098896194e-07,
        -7.37767982771374e-07,
        -9.90342442314471e-07,
        -1.25224526326335e-06,
        -1.50325781867309e-06,
        -1.71932096193714e-06,
        -1.87278144550817e-06,
        -1.93262925122866e-06,
        -1.86472603152310e-06,
        -1.63202495673730e-06,
        -1.19478222735111e-06,
        -5.10760494005637e-07]

        @test all(hrf_ref .≈ hrf[1:20])

        dcm.U = nothing
        hrf = rDCM.get_hrf(dcm)
        hrf_ref = [0.0,
            0.0,
            0.0,
            0.0,
            -3.291034022872453e-7,
            -1.1186551508966146e-6,
            -2.323319830470751e-6,
            -3.740689203218181e-6,
            -5.029503350365019e-6,
            -5.726565529015204e-6,
            -5.26242693096037e-6,
            -2.975913690464425e-6,
            1.8724351421153048e-6,
            1.0087953622795571e-5,
            2.2530397133001847e-5,
            4.0103394435775375e-5,
            6.37447482576026e-5,
            9.441753206153054e-5,
            0.00013310193361654673,
            0.00018078779859944974]
        @test all(hrf_ref .≈ hrf[1:20])

        dcm.Y.y = nothing
        @test_throws ErrorException("BOLD signal is empty.") rDCM.get_hrf(dcm)
    end
end

function test_generate_linear(dcm)

    @testset "linear" begin
        _, y, x, _ = generate_BOLD(dcm;SNR=3)

        y_ref = [-0.0495902545634915,
        -0.0429054182953994,
        -0.0370705932010462,
        -0.0319854078525914,
        -0.0275602832310189,
        -0.0237154401097601,
        -0.0203799883352225,
        -0.0174911530915810,
        -0.0149937006836198,
        -0.0128396105383830,
        -0.0109879980445423,
        -0.00940523459418523,
        -0.00806515374812553,
        -0.00694919238276097,
        -0.00604630393582438,
        -0.00535250062304212,
        -0.00486992864358030,
        -0.00460544567003442,
        -0.00456874164626957,
        -0.00477011055803912]

        x_ref = [-0.000355676052658233,
        -0.000352092210928230,
        -0.000348543162206672,
        -0.000345028584953782,
        -0.000341548160401681,
        -0.000338101572532960,
        -0.000334688508059374,
        -0.000331308656400682,
        -0.000327961709699946,
        -0.000324647362872554,
        -0.000321365313686301,
        -0.000318115262870018,
        -0.000314896914248332,
        -0.000311709974900227,
        -0.000308554155339226,
        -0.000305429169713059,
        -0.000302334736020817,
        -0.000299270576345669,
        -0.000296236417101303,
        -0.000293231989290359]

        @test all(y_ref .≈ y[1:20,1])
        @test all(x_ref .≈ x[1:20,1])
    end
    @testset "special cases" begin
        @test_logs (:warn,"Overwriting sampling time of Y with TR.") generate_BOLD(dcm;SNR=3,TR=1.0)
        dcm.Y = nothing
        generate_BOLD(dcm;SNR=3, TR=1.0) #TODO: make test for results
        dcm = load_example_DCM()
        dcm.U.u .*= 10.0 # set u to larger values such that BOLD signal values go above 20
        @test_logs (:warn,"BOLD signal contains large values.") generate_BOLD(dcm;SNR=3)
        dcm.Ep.A[1,1] = 10 # make system unstable
        @test_logs (:warn,"BOLD signal contains Inf or NaN. Check data generating parameters or input structure.") generate_BOLD(dcm;SNR=3)
        dcm.U.u[:] .= 0.0 # no input that would drive the system
        @test_logs (:warn,"BOLD signal contains only zeros. Check data generating parameters or input structure.") generate_BOLD(dcm;SNR=3)
    end
end

function test_generate_bilinear()
    @testset "bilinear" begin
        y_ref = [0.0321329949386491,
        0.0663211733707552,
        -0.0496818678617098,
        -0.0291297937814116,
        0.0491892650361192,
        0.000547377580766762,
        -0.146991221824794,
        -0.470753548170548,
        -0.596179680668909,
        -0.360149042362737,
        -0.204353198143214,
        -0.287256262868103,
        -0.393299926855943,
        -0.280637531632456,
        0.0679882602052085,
        0.175466106087251,
        -0.0250102972310988,
        -0.0644816067067507,
        0.0501738458750545,
        0.0161010418455187]

        x_ref = [-0.0154453641179358,
        	-0.0144256000229186,
            -0.0132396170894261,
            -0.0119054496884601,
            -0.0104427075772449,
            -0.00887229820178303,
            -0.00721613729585289,
            -0.00549685167298369,
            -0.00373747813946321,
            -0.00196116243729003,
            -0.000190862057305801,
            0.00155094335372603,
            0.00324253043026921,
            0.00486319601211661,
            0.00639350344197849,
            0.00781550840235930,
            0.00911296176807220,
            0.0102714872880572,
            0.0112787322694846,
            0.0121244898134723]

        dcm_bi = rDCM.load_DCM_bilinear()
        _, y, x, _ = generate_BOLD(dcm_bi;SNR=3)

        @test all(y_ref .≈ y[1:20,1])
        @test all(x_ref .≈ x[1:20,1])
    end
end

function test_generate_nonlinear()
    @testset "nonlinear" begin
        y_ref = [0.0330336291743252,
        0.0665343424300450,
        -0.0456961116174254,
        -0.0272865307283921,
        0.0481491930882374,
        0.00164387155647311,
        -0.144927339313638,
        -0.464993490279593,
        -0.574813789380730,
        -0.341123821613434,
        -0.206535004740480,
        -0.291082929712693,
        -0.381652114399425,
        -0.267514459311530,
        0.0700034424041667,
        0.178888565176907,
        -0.0175121522928192,
        -0.0628681644556771,
        0.0488851979895122,
        0.0192157183464367]

        x_ref = [-0.0149604089047334,
        -0.0139847279152765,
        -0.0128437849121769,
        -0.0115556998472636,
        -0.0101402510696320,
        -0.00861855109053162,
        -0.00701270987706510,
        -0.00534549254639817,
        -0.00363997810210866,
        -0.00191922543555109,
        -0.000205952239704804,
        0.00147776821260329,
        0.00311078825149422,
        0.00467314160008239,
        0.00614627124653424,
        0.00751322767712344,
        0.00875883650215976,
        0.00986983518689852,
        0.0108349791917235,
        0.0116451183205517]

        dcm_non = rDCM.load_DCM_nonlinear()
        _, y, x, _ = generate_BOLD(dcm_non;SNR=3)

        @test all(y_ref .≈ y[1:20,1])
        @test all(x_ref .≈ x[1:20,1])
    end
end

function test_generate_pureJulia()

    y_ref = [0.0330336291743252,
    0.0665343424300450,
    -0.0456961116174254,
    -0.0272865307283921,
    0.0481491930882374,
    0.00164387155647311,
    -0.144927339313638,
    -0.464993490279593,
    -0.574813789380730,
    -0.341123821613434,
    -0.206535004740480,
    -0.291082929712693,
    -0.381652114399425,
    -0.267514459311530,
    0.0700034424041667,
    0.178888565176907,
    -0.0175121522928192,
    -0.0628681644556771,
    0.0488851979895122,
    0.0192157183464367]

    x_ref = [-0.0149604089047334,
    -0.0139847279152765,
    -0.0128437849121769,
    -0.0115556998472636,
    -0.0101402510696320,
    -0.00861855109053162,
    -0.00701270987706510,
    -0.00534549254639817,
    -0.00363997810210866,
    -0.00191922543555109,
    -0.000205952239704804,
    0.00147776821260329,
    0.00311078825149422,
    0.00467314160008239,
    0.00614627124653424,
    0.00751322767712344,
    0.00875883650215976,
    0.00986983518689852,
    0.0108349791917235,
    0.0116451183205517]

    # rename binary file such that package is forced to use Julia implementation
    orig_name = joinpath(rDCM.euler_integration_bin, "libdcm_euler_integration.so")
    alt_name = joinpath(rDCM.euler_integration_bin, "hidden.so")
    mv(orig_name,alt_name)

    dcm_non = rDCM.load_DCM_nonlinear()
    _, y, x, _ = generate_BOLD(dcm_non;SNR=3)

    #rename file to original name
    mv(alt_name,orig_name)
    @test all(y_ref .≈ y[1:20,1])
    @test all(x_ref .≈ x[1:20,1])

end

function test_error_handling(dcm)

    # specified TR is not a multiple of u_dt
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
            of the input U (u_dt). Cannot proceed.") generate_BOLD(dcm;SNR=3,TR=0.66)

    # sampling rate of Y is not a multiple of sampling rate of U
    dcm.Y.dt = 0.66
    @test_throws ErrorException("The sampling rate of Y (y_dt) is not a multiple of the sampling rate
                  of the input U (u_dt). Cannot proceed.") generate_BOLD(dcm;SNR=3)

    # no sampling rate provided
    dcm.Y = nothing
    @test_throws ErrorException("Y field is empty. Please specify the TR of the BOLD signal manually
            using tapas_rdcm_generate(dcm;TR=value)") generate_BOLD(dcm;SNR=3)

end

function test_simulate_data()

    dcm = load_example_DCM()

    test_HRF(copy(dcm))
    @testset verbose=true "Generate synthetic data" begin
        test_generate_linear(copy(dcm))
        test_error_handling(copy(dcm))
        test_generate_bilinear()
        test_generate_nonlinear()
        test_generate_pureJulia()
    end
end

test_simulate_data()
