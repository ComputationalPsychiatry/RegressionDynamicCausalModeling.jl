
function predict_rigid_rdcm(dcm)
    @testset "rigid rDCM" begin

        # set options for inversion
        opt = Options(RigidInversionParams();synthetic=true,
        verbose=0,
        testing=true)
        rdcm = RigidRdcm(dcm)

        output = invert(rdcm, opt)

        y_pred = predict(rdcm,output)

        # freq domain result
        # y_pred_ref = [
        #     -0.874692922856965,
        #     -0.877473052415406,
        #     -0.8963963522140176,
        #     -0.9318195308449948,
        #     -0.983778139932242,
        #     -1.0519902020764298,
        #     -1.1358722588349424,
        #     -1.234565832181771,
        #     -1.346971975216602,
        #     -1.471791482122119
        # ]

        # time domain result
        y_pred_ref = [
            -0.8450826325060272,
            -0.8492376534030887,
            -0.8692779017868778,
            -0.9055312501368019,
            -0.9580081075449519,
            -1.0264058275139354,
            -1.1101253808057039,
            -1.2082982610439015,
            -1.319821285996252,
            -1.4433968635718357
        ]

        y_idx = 1300

        @test all(y_pred_ref .≈ y_pred[y_idx+1:y_idx+10,1])
    end
end

function predict_sparse_rdcm(dcm)
    @testset "sparse rDCM" begin
        # create rDCM struct
        rdcm = SparseRdcm(dcm;p0=0.15)

        # set options for inversion
        opt = Options(SparseInversionParams(;reruns=10,restrictInputs=true);
            synthetic=true,
            verbose=0,
            testing=true)

        output = invert(rdcm,opt)

        y_pred = predict(rdcm,output)

        # freq domain result
        # y_pred_ref = [
        #     -0.8149854390907948,
        #     -0.8011610787824592,
        #     -0.8017856907878724,
        #     -0.8174094279497025,
        #     -0.8483001943993492,
        #     -0.8944332007979571,
        #     -0.9554935449488213,
        #     -1.0308908634969143,
        #     -1.1197843839436636,
        #     -1.2211162138586589
        # ]

        # time domain result
        y_pred_ref = [
            -0.9163667995842838,
            -0.9216556567351829,
            -0.9425146052457634,
            -0.9791877648486933,
            -1.0316056351344263,
            -1.0993955786683156,
            -1.181903853404036,
            -1.2782267654592425,
            -1.3872483478678481,
            -1.5076820182918629
        ]

        y_idx = 1300

        @test all(y_pred_ref .≈ y_pred[y_idx+1:y_idx+10,1])
    end
end

function predict_error_handling(dcm)

    @testset "error handling" begin

        # set options for inversion
        opt = Options(RigidInversionParams();synthetic=true,
        verbose=0,
        testing=true)
        rdcm = RigidRdcm(dcm)

        output = invert(rdcm, opt)
        rdcm.c = BitMatrix(zeros(size(rdcm.c)))

        @test_throws ErrorException("Cannot generate data from resting-state DCM.") predict(rdcm,output)

    end
end

function test_prediction()
    dcm = load_example_DCM()
    @testset verbose=true "Prediction" begin
        predict_rigid_rdcm(copy(dcm))
        predict_sparse_rdcm(copy(dcm))
        predict_error_handling(copy(dcm))
    end
end

test_prediction()
