
function test_rigidInversion()

    @testset "rigid rDCM" begin

        dcm = load_example_DCM()

        # create rDCM struct
        dcm = RigidRdcm(dcm)

        # set options for inversion
        opt = Options(RigidInversionParams();synthetic=true,verbose=0,testing=true)

        output = invert(dcm,opt)

        # F_ref = 34658.2168411653 # result from Matlab
        F_ref = 75495.9676969727 # Mathematically correct result
        a_ref = zeros(50) .+ 1359.0

        # Matlab result
        # b_ref = [278.992934607013,
        # 279.930331980155,
        # 1232.80824967765,
        # 68.7784197407701,
        # 383.400375214554,
        # 196.978150642330,
        # 504.246840658795,
        # 212.196814776145,
        # 51.1264041635164,
        # 72.4320581916030,
        # 77.9692348000759,
        # 78.4644701413868,
        # 263.955567818095,
        # 11.8783419613710,
        # 20.7426505521006,
        # 87.4782008875842,
        # 401.376875928776,
        # 569.781344435128,
        # 4.78095908170428,
        # 12.4220174434528]

        # m_ref = [-0.416024957356795	0.395111621322398	0	0	0.288189637587704	0	0	0	0;
        # 0	-0.524526786993822	-0.271165356082288	0	0	0	0	0	0;
        # 0	0	-0.465482610530771	0.750847061568657	0	0	0	0.797990342079961	0;
        # 0	0	0	-0.470684309752512	0.257011870842332	0	0	0	0;
        # 0	0	0	0	-0.446990648378351	0	0	0	0;
        # 0	0	0	0	0	-0.523978201475677	0.363189573059244	0	0;
        # 0	0	0	0	0	0	-0.542823940020913	0.906099594873587	0;
        # 0	0	0	0	0	0	0	-0.506244060272773	-0.560825878389055;
        # 0	0	0	0	0	0	0	0	-0.538796358363997]

        # Σ_ref = [1.22710745970774e-07	-1.10981038945065e-07	0	0	-8.67497517499543e-08;
        # -1.10981038945065e-07	1.07312798959610e-07	0	0	8.17805166976585e-08;
        # 0	0	0	0	0;
        # 0	0	0	0	0;
        # -8.67497517499542e-08	8.17805166976585e-08	0	0	7.68038324402694e-08]

        # Mathematically correct results
        b_ref= [144.4274213518344,
        154.7393100385269,
        679.9853727258028,
        36.22214576430868,
        198.1069902794391,
        102.49195105002349,
        266.9653607889319,
        114.8333495760788,
        27.031062340259364,
        37.83635099545427,
        39.91228451387678,
        42.577718804528516,
        140.28033756436812,
        6.647262969053058,
        11.175953980319779,
        46.826165270045145,
        217.33288740410183,
        302.5076422861302,
        2.958767451047271,
        6.879805374420489]

        m_ref = [-0.41548922171677133 0.3949738877798173 0.0 0.0 0.287735406093434 0.0 0.0 0.0 0.0;
        0.0 -0.5238266884300515 -0.2711682814395596 0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 -0.46454372168894786 0.750890653792754 0.0 0.0 0.0 0.7973417913312492 0.0;
        0.0 0.0 0.0 -0.47034688567970767 0.25701330477855067 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 -0.4466242012603222 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 -0.5233253247330675 0.3629654752537566 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0 -0.5422371148997025 0.9058915113714406 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.5055170998094115 -0.5606792333188154;
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.5384563384361282]

        Σ_ref = [1.2359990303755105e-7 -1.1401652455841153e-7 0.0 0.0 -8.689259265429028e-8;
        -1.1401652455841153e-7 1.1087742770185532e-7 0.0 0.0 8.39189348225854e-8;
        0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0;
        -8.689259265429027e-8 8.391893482258539e-8 0.0 0.0 7.704000346831523e-8]

        @test F_ref ≈ output.F
        @test all(a_ref .≈ output.a_all)
        @test all(b_ref .≈ output.b_all[1:20])
        @test all(m_ref .≈ output.m_all[1:9,1:9])
        @test all(Σ_ref .≈ output.Σ_all[1][1:5,1:5])

        # Matlab result
        # A_mat = [-0.416024957356795	0.395111621322398	0	0	0.288189637587704;
        # 0	-0.524526786993822	-0.271165356082288	0	0;
        # 0	0	-0.465482610530771	0.750847061568657	0;
        # 0	0	0	-0.470684309752512	0.257011870842332;
        # 0	0	0	0	-0.446990648378351]

        # C_mat = [0	-0.568961980009389	0	0	0	0;
        # 0	0	0	0	0	0;
        # 0	0	0.942108612362386	0	0	0;
        # 0	0	0	0	0	0;
        # 0	0	0	0	1.96255418421427	0;
        # 0	0	0	0	0	0;
        # 0	0	0	0	0	-0.579412563427436]

        A_mat = [-0.41548922171677133 0.3949738877798173 0.0 0.0 0.287735406093434;
        0.0 -0.5238266884300515 -0.2711682814395596 0.0 0.0;
        0.0 0.0 -0.46454372168894786 0.750890653792754 0.0;
        0.0 0.0 0.0 -0.47034688567970767 0.25701330477855067;
        0.0 0.0 0.0 0.0 -0.4466242012603222]

        C_mat = [0.0 -0.5690155750251027 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.9427587266691153 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 1.9625660471314932 0.0;
        0.0 0.0 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 -0.5771997855126901]

        out_dcm = LinearDCM(dcm,output)

        @test all(A_mat .≈ out_dcm.Ep.A[1:5,1:5])
        @test all(C_mat .≈ out_dcm.Ep.C[1:7,1:6])
    end
end

function test_sparseInversion()

    @testset "sparse rDCM (fix C matrix)" begin

        dcm = load_example_DCM()

        # create rDCM struct
        dcm = SparseRdcm(dcm;p0=0.15)

        # set options for inversion
        opt = Options(SparseInversionParams(;reruns=10,restrictInputs=true);
            synthetic=true,
            verbose=0,
            testing=true)

        out = invert(dcm,opt)

        # original values from matlab. Region 7 differs (thats where negative b occurs in matlab implementation)
        # F_ref = [-683.015071423908;
        # 	81.3162140763429;
        #     -7040.98079782687;
        #     980.452205220891;
        #     -152.404136020184;
        #     -7539.81777857265;
        #     1023.07951193787;
        #     2508.10599848526;
        #     1512.64366485487;
        #     856.177638597529]

        # Matlab version but taking the real part after every update equation
        # F_ref = [-683.015071432846;
        # 	81.3162140735612;
        #     -7040.98079782688;
        #     980.452205223091;
        #     -152.404136023463;
        #     -7539.81777857266;
        #     1032.80577292814;
        #     2508.10599849728;
        #     1512.64366483464;
        #     856.177638596816]

        a_ref = zeros(50) .+ 1359.0

        # original values from matlab.
        # b_ref = [201.856976601018,
        # 113.731780446923,
        # 22039.7335136045,
        # 56.2417795222879,
        # 144.217765295817,
        # 31422.8145016882,
        # 55.8163377536835,
        # 18.2416329222937,
        # 40.4877739115726,
        # 68.2924505508023,
        # 3679.56879215311,
        # 1994.10097846738,
        # 139.969172314735,
        # 11.2957605003429,
        # 20.7426505526101,
        # 2106.32796178814,
        # 196.388889047223,
        # 13794.5029572736,
        # 4.72511578237763,
        # 12.4220174438525]

        # Matlab version but taking the real part after every update equation
        # b_ref = [201.856976602262,
        # 113.731780446563,
        # 22039.7335136052,
        # 56.2417795223178,
        # 144.217765296167,
        # 31422.8145016880,
        # 55.8447556415517,
        # 18.2416329237142,
        # 40.4877739115617,
        # 68.2924505508387,
        # 3679.56879215313,
        # 1.11473277645656,
        # 139.969172315345,
        # 11.2957605003620,
        # 20.7426505525965,
        # 2106.32796178795,
        # 196.388889048544,
        # 13794.5029572735,
        # 4.72511578238127,
        # 12.4220174438511]

        # original values from matlab.
        # m_ref = [-0.410165273077682	0.387946381496301	0	0	0.283867214559243	0	0.00213784031872032	0	0.00711292709137987;
        # 0.00285962036408516	-0.527629151460803	-0.268563178105413	-0.00797591531107030	0	-0.0256275481484689	0.0211489952794854	0.00170052149143193	0.0396055338420940;
        # -0.00581861173924450	0	-0.493910058573326	0.801012561347083	0	0	0	0.818407380376522	-0.0399525596455128;
        # -0.000728337151827026	0	0	-0.471586378087789	0.257728588557155	0	0	0	0.00205091386540074;
        # 0	0	0	0.110709366643499	-0.550501799656266	0	0	0	0;
        # 0.0161175008163074	0.151541702482121	0.0759229797007047	0	0	-0.0620329236998051	0	0.0923887706928968	0;
        # 0	0	7.49607874108281e-05	0	0	0.0946347383355348	-0.659538277295524	0.989694849859056	-0.0917656062422947;
        # -0.000266752301711002	0	0	0	0	0.00233043077765599	0	-0.477730394001841	-0.486597929486026;
        # 0.000177524276320460	0	0	0	0	-0.00251421792880777	0	0	-0.545619327965552]

        # Matlab version but taking the real part after every update equation
        # m_ref = [-0.410165273077681	0.387946381496298	0	0	0.283867214559242	0	0.00213784031871940	0	0.00711292709137967;
        # 0.00285962036408611	-0.527629151460815	-0.268563178105423	-0.00797591531106226	0	-0.0256275481484664	0.0211489952794768	0.00170052149145459	0.0396055338420905;
        # -0.00581861173925369	0	-0.493910058573355	0.801012561347129	0	0	0	0.818407380376584	-0.0399525596454816;
        # -0.000728337151827407	0	0	-0.471586378087787	0.257728588557153	0	0	0	0.00205091386540074;
        # 0	0	0	0.110709366643524	-0.550501799656292	0	0	0	0;
        # 0.0161175008163080	0.151541702482119	0.0759229797007043	0	0	-0.0620329236998058	0	0.0923887706928972	0;
        # 0	0	0	0	0	0.0947317419078151	-0.659620944485194	0.989938866346546	-0.0917685355815225;
        # -0.000266752301710606	0	0	0	0	0.00233043077766015	0	-0.477730394001949	-0.486597929486148;
        # 0.000177524276320718	0	0	0	0	-0.00251421792880737	0	0	-0.545619327965550]

        # Matlab version but taking the real part after every update equation
        # Σ_ref = [1.33840097727108e-07	-1.30849312093877e-07	0	0	-9.75251190312053e-08;
        # -1.30849312093877e-07	1.79405305731737e-07	0	0	1.11823736037588e-07;
        # 0	0	0.160000000000000	0	0;
        # 0	0	0	0.160000000000000	0;
        # -9.75251190312053e-08	1.11823736037588e-07	0	0	8.95973250213253e-08]

        # original values from matlab.
        # z_ref = [1	1	0;
        # 1	1	1;
        # 1.00000000000000	0	1;
        # 1	0	0;
        # 0	0	0;
        # 1	1	1;
        # 0	0	0.999956312318118;
        # 1	0	0;
        # 0.999999672922171	0	0;
        # 1	0	0;
        # 0.999999999082709	1	0]

        # Matlab version but taking the real part after every update equation
        # z_ref = [1	1	0;
        # 1	1	1;
        # 1.00000000000000	0	1;
        # 1	0	0;
        # 0	0	0;
        # 1	1	1;
        # 0	0	0;
        # 1	0	0;
        # 0.999999672922171	0	0;
        # 1	0	0;
        # 0.999999999082709	1	0]

        # Result if QF multiplied extra with 0.5 as in Matlab (which is wrong)
        # F_ref = [227.51060418136126,
        # 817.2585668200674,
        # -6108.053655687886,
        # 1922.4353674234837,
        # 775.3172570485297,
        # -6587.381679658058,
        # 1931.796078956848,
        # 3349.703873466042,
        # 2400.441347901507,
        # 1765.1581578681094]

        F_ref = [-445.39116404238706,
        147.92914454706982,
        -6787.4919688939735,
        1266.5519849880081,
        105.17446118685854,
        -7266.838113628709,
        1275.9012095972055,
        2740.6985525309756,
        1753.5141906834,
        1105.220119583982]

        # Result if QF multiplied extra with 0.5 as in Matlab (which is wrong)
        # b_ref = [102.98213573071136,
        # 66.81036399691982,
        # 11015.3236520329,
        # 28.772113455430333,
        # 72.61784502666174,
        # 15597.014695847918,
        # 28.786114779566944,
        # 9.639025361774967,
        # 20.86093634964251,
        # 34.73578007509626,
        # 1845.6205942712604,
        # 1.0579624264226333,
        # 12964.540406017166,
        # 6.377854819974421,
        # 11.175953980766646,
        # 84.11547718288665,
        # 103.01108580008145,
        # 6196.216671355365,
        # 2.8860123258329633,
        # 6.879805374281059]

        b_ref = [102.98213812266839,
        66.81036399691982,
        11015.3236520329,
        28.772113455430333,
        72.61784502666174,
        15597.014695847918,
        28.786114779566944,
        9.639025380095969,
        20.86093634964251,
        34.73578007509626,
        1845.6205942712604,
        1.05796242684734,
        12964.540406017166,
        6.377854819974421,
        11.175953980766646,
        84.11547718288665,
        103.01108580008145,
        6196.216671355365,
        2.8860123258329633,
        6.879805374281059]

        m_ref = [-0.4105086155216504 0.3889560700865953 0.0 0.0 0.2843449533998572 0.0 0.002634915014976799 0.0 0.007515301587383255;
        0.003984618745544368 -0.5314141001437789 -0.27050074548890624 -0.008670121696506597 0.0 -0.028266541494481073 0.022717749900146595 0.004403665037303528 0.04383754373674989;
        0.0028757328887315525 0.0 -0.4690265232823516 0.6963871817766307 0.033969160213959015 0.0 0.0 0.7149956488927196 -0.11456351237630247;
        -0.0007566714171183403 0.0 0.0 -0.4712446035286818 0.25761104444588345 0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.11090542021157734 -0.5506713830892727 0.0 0.0 0.0 0.0; 0.014357127711714174 0.18436211941261899 0.09473980035927192 0.0 0.0 -0.057971623014369686 0.0 0.08740921255381469 0.0;
        0.0 0.0 6.795149088469383e-5 0.0 0.0 0.09571957209998172 -0.6608752128148976 0.9909784437500648 -0.09245981846928952;
        -0.0002769500071787105 0.0 0.0 0.0 0.0 0.00240276167681458 0.0 -0.4775922647774866 -0.48614355299834866;
        0.0002795873572726131 0.0 0.0 0.0 0.0 -0.0024627625295305894 0.0 0.0 -0.5456422803429525]

        # Result if QF multiplied extra with 0.5 as in Matlab (which is wrong)
        # Σ_ref = [1.3238005689503363e-7 -1.3023688352337991e-7 0.0 0.0 -9.600712932694173e-8;
        # -1.3023688352338018e-7 1.7615197467920715e-7 0.0 0.0 1.1013550032456448e-7;
        # 0.0 0.0 0.16 0.0 0.0; 0.0 0.0 0.0 0.16 0.0;
        # -9.600712932694175e-8 1.1013550032456435e-7 0.0 0.0 8.824935711952757e-8]

        Σ_ref = [1.323806141079705e-7 -1.302374317110312e-7 0.0 0.0 -9.600753343754925e-8;
        -1.3023743171103146e-7 1.761527161372645e-7 0.0 0.0 1.1013596390386595e-7;
        0.0 0.0 0.16 0.0 0.0; 0.0 0.0 0.0 0.16 0.0;
        -9.600753343754934e-8 1.1013596390386587e-7 0.0 0.0 8.82497285800257e-8]

        z_ref = [1.0 1.0 0.0;
        1.0 1.0 1.0;
        0.9999999962943193 0.0 1.0;
        1.0 0.0 0.0;
        0.0 0.0 0.0;
        1.0 1.0 1.0;
        0.0 0.0 0.9999675299265609;
        1.0 0.0 0.0;
        1.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0]

        @test all(F_ref .≈ out.F_r[1:10])
        @test all(a_ref .≈ out.a_all)
        @test all(b_ref .≈ out.b_all[1:20])
        @test all(m_ref .≈ out.m_all[1:9,1:9])
        @test all(Σ_ref .≈ out.Σ_all[1][1:5,1:5])
        @test all(z_ref .≈ out.z_all[1:11,1:3])


    end

    @testset "sparse rDCM (prune C matrix)" begin

        dcm = load_example_DCM()

        # create rDCM struct
        dcm = SparseRdcm(dcm;p0=0.15)

        # set options for inversion
        opt = Options(SparseInversionParams(;reruns=10,restrictInputs=false);
            synthetic=true,
            verbose=0,
            testing=true)

        out = invert(dcm,opt)

        F_ref = [603.1616658254094,
            853.1585281139579,
            -1094.348026036079,
            5674.542786572273,
            -15.38888788470415,
            -6916.189922604427,
            -5802.550192324399,
            -3640.7540663495074,
            5834.670014640793,
            2397.1028420119933]

        a_ref = zeros(50) .+ 1359.0

        b_ref = [44.647805284377014,
            38.31464306370982,
            158.95186687483599,
            1.0944952081179433,
            72.39768053695335,
            11803.740437559167,
            4790.440786524036,
            1027.7603427014692,
            1.0693040041383075,
            12.781594080832114,
            36.3894888778394,
            6.51761377921766,
            179.89759196847785,
            5.9765977100376375,
            4.4322394680157515,
            4538.957587205635,
            1.362796355440854,
            4848.551523707802,
            1.0069470406052246,
            2.8087200713262668]

        m_ref = [-0.4033266445569739 0.3811485384970395 0.0 0.055370112802457974;
            0.016004087591951377 -0.5492792077069786 -0.27441471831652753 -0.02336387803955281;
            -0.03953284781938693 0.0 -0.5107470886146426 0.8450109730721679;
            0.0 0.0 0.0 -0.4157776986122256]

        Σ_ref = [1.3811302976845223e-7 -1.600638052793148e-7 0.0 -1.5948331476051038e-7;
            -1.6006380527931415e-7 2.2600468569930597e-7 0.0 3.2214512088137197e-7;
            0.0 0.0 0.16 0.0;
            -1.5948331476049816e-7 3.2214512088135233e-7 0.0 2.6456519301683646e-6]

        z_ref = [1.0 1.0 0.0;
            1.0 1.0 1.0;
            1.0 0.0 1.0;
            0.0 0.0 0.0;
            1.0 1.0 0.0;
            1.0 1.0 1.0;
            1.0 0.0 0.9999992503923544;
            1.0 0.0 0.0;
            0.0 0.0 0.0;
            0.9981699204944517 0.0 0.0;
            0.9999999999999905 1.0 1.0]

        @test all(F_ref .≈ out.F_r[1:10])
        @test all(a_ref .≈ out.a_all)
        @test all(b_ref .≈ out.b_all[1:20])
        @test all(m_ref .≈ out.m_all[1:4,1:4])
        @test all(Σ_ref .≈ out.Σ_all[1][1:4,1:4])
        @test all(z_ref .≈ out.z_all[1:11,1:3])
    end

end

function test_model_inversion()
    @testset verbose=true "Model inversion" begin
        test_rigidInversion()
        test_sparseInversion()
    end
end

test_model_inversion()
