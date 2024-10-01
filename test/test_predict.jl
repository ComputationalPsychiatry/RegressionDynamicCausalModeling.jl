
function predict_rigid_rdcm(dcm)
    @testset "rigid rDCM" begin

        # set options for inversion
        opt = Options(RigidInversionParams();synthetic=true,
        verbose=0,
        testing=true)
        rdcm = RigidRdcm(dcm)

        output = invert(rdcm, opt)

        y_pred = predict(rdcm,output)

        # prediction of Matlab implementation
        # y_pred_ref = [-0.0462946041496438	-0.0149355589977407;
        # -0.0400225946027140	-0.0126340112965214;
        # -0.0345525116007612	-0.0106681743193625;
        # -0.0297890241474584	-0.00899245583814413;
        # -0.0256471531126566	-0.00756711780022254;
        # -0.0220513012368967	-0.00635757750760177;
        # -0.0189343667176012	-0.00533402638443777;
        # -0.0162369929582907	-0.00447149206754289;
        # -0.0139070142348885	-0.00375039026458674;
        # -0.0118991415370625	-0.00315746228475652;
        # -0.0101748919168788	-0.00268684039066030;
        # -0.00870270823757472	-0.00234088738236362;
        # -0.00745816081773444	-0.00213045005222443;
        # -0.00642408416767119	-0.00207424667609078;
        # -0.00559049132223353	-0.00219725154426424;
        # -0.00495412811020659	-0.00252810901910400;
        # -0.00451757592662251	-0.00309577042865277;
        # -0.00428787509285637	-0.00392567166553826;
        # -0.00427471031626602	-0.00503584120679113;
        # -0.00448826392413132	-0.00643334239519629]

        y_pred_ref = [-0.04701485434644018 -0.01516880174832442;
        -0.040655513996814954 -0.012834709150314559;
        -0.03510782832192189 -0.010840622457606574;
        -0.03027551793080133 -0.009140418398939596;
        -0.026072725142971185 -0.007693888998987735;
        -0.0224230432887253 -0.0064660365982712985;
        -0.019258628502807573 -0.005426687423197247;
        -0.01651944702177433 -0.004550550004246451;
        -0.014152718190900398 -0.003817765340872802;
        -0.012112597872095641 -0.003214844797958642;
        -0.010360105923745613 -0.0027357374326079345;
        -0.008863244806708533 -0.0023826721981328775;
        -0.007597200739770834 -0.002166413358893994;
        -0.006544480287005988 -0.0021056479023473127;
        -0.005694824368956448 -0.0022253667807017766;
        -0.0050447614144659755 -0.0025542717388085086;
        -0.004596707596574964 -0.003121400888146985;
        -0.004357585727889246 -0.003952291359360457;
        -0.004337004025134855 -0.00506506975696424;
        -0.004545100383626692 -0.006466875665109336]

        @test all(y_pred_ref .≈ y_pred[1:20,1:2])
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

        # Result if QF multiplied extra with 0.5 as in Matlab (which is wrong)
        # y_pred_ref = [-0.024646337704718926 -0.0013352067108216175;
        # -0.02020948152024477 3.790719868965559e-6;
        # -0.016380879837471447 0.0011464263653869317;
        # -0.013080659238681057 0.0021158016098709056;
        # -0.010243543008730776 0.0029312422447245346;
        # -0.007820935010489072 0.0036098640056828517;
        # -0.005781357258266652 0.004167552415960477;
        # -0.004108064608630804 0.00461857837159202;
        # -0.0027944180758166116 0.004973520212353151;
        # -0.0018384678399396607 0.005235861145427054;
        # -0.0012381977286864776 0.005398177607577026;
        # -0.0009883852482597758 0.00543903797381928;
        # -0.0010793936945394055 0.005321591054925603;
        # -0.0014976701280489254 0.00499444773020703;
        # -0.0022273819151109177 0.004394986571539735;
        # -0.00325249870399105 0.003454767810519021;
        # -0.004558676193453195 0.002106401899673764;
        # -0.006134461406473561 0.00029102814932705065;
        # -0.007971554521291838 -0.00203448192969999;
        # -0.010064078824555248 -0.004891395301152602]

        y_pred_ref = [-0.02464636623847387 -0.0013352258696505586;
        -0.020209507949086866 3.7733382105444027e-6;
        -0.016380904250388684 0.001146410650061764;
        -0.01308068172750108 0.0021157874520108446;
        -0.010243563667205908 0.002931229538054421;
        -0.007820953933508948 0.003609852646941079;
        -0.005781374541282905 0.004167542304736349;
        -0.004108080347382722 0.0046185694089158815;
        -0.0027944323665640746 0.004973512297198092;
        -0.0018384807804727313 0.0052358541684647064;
        -0.0012382094204109384 0.005398171442007971;
        -0.000988395799611786 0.005439032463751804;
        -0.0010794032257957574 0.005321586002421454;
        -0.0014976787773969618 0.004994442882555644;
        -0.002227389845572565 0.00439498161053612;
        -0.0032525061105865135 0.0034547623453954623;
        -0.00455868330984995 0.0021063954649512746;
        -0.006134468510318826 0.0002910202083456148;
        -0.007971561937487905 -0.0020344919754783104;
        -0.010064086925918233 -0.004891408097538302]

        @test all(y_pred_ref .≈ y_pred[1:20,1:2])
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
