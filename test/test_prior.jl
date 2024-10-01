
function test_pE(a::BitArray,c::BitArray)

  @testset "Prior param" begin

    A_pE = [-0.5   0.0   0.0
    0.0  -0.5   0.0
    0.0   0.0  -0.5]

    C_pE = [0	0	0	0;
    0	0	0	0;
    0	0	0	0]

    A_pC = [0.0416666666666667	0	2.66666666666667;
    2.66666666666667	0.0416666666666667	0;
    0	0	0.0416666666666667]

    C_pC = [1	1	1	1;
    0	0	1	1;
    1	0	1	0]

    transit_pE = zeros(3)
    decay_pE = zeros(3)
    epsilon_pE = 0

    transit_pC = [0.00247875217666636,
    0.00247875217666636,
    0.00247875217666636]
    decay_pC = [0.00247875217666636,
    0.00247875217666636,
    0.00247875217666636]
    epsilon_pC = 0.00247875217666636

    prior = rDCM.get_prior_stats(a,c)

    @test all(A_pE .≈ prior.pE.A)
    @test all(C_pE .≈ prior.pE.C)
    @test all(transit_pE .≈ prior.pE.transit)
    @test all(decay_pE .≈ prior.pE.decay)
    @test all(epsilon_pE .≈ prior.pE.epsilon)

    @test all(A_pC .≈ prior.pC.A)
    @test all(C_pC .≈ prior.pC.C)
    @test all(transit_pC .≈ prior.pC.transit)
    @test all(decay_pC .≈ prior.pC.decay)
    @test all(epsilon_pC .≈ prior.pC.epsilon)
  end
end

function test_sample_prior(a::BitArray,c::BitArray)

  @testset "Sample prior" begin
    # construct a test where the system is unstable and check if function detects that
    prior = rDCM.get_prior_stats(a,c)
    # create a prior mean such that system will be unstable
    A_unstable = [0.5 0.0 0.0052;
    0.0052 0.5 0.0;
    0.0 0.0 0.5]

    pE = rDCM.priorMeanLinear(A_unstable,prior.pE.C,prior.pE.transit,prior.pE.decay,prior.pE.epsilon)
    pC = rDCM.priorCovLinear(prior.pC.A,prior.pC.C,prior.pC.transit,prior.pC.decay,prior.pC.epsilon)

    @test_throws ErrorException("Not able so sample values such that system is stable.") sample_from_prior(a,c,rDCM.PriorDCMLinear(pE,pC);rng=MersenneTwister(rDCM.FIXEDSEED))

    # normal example
    prior = rDCM.get_prior_stats(a,c)
    A, C, transit, decay, epsilon = sample_from_prior(a,c,prior;rng=MersenneTwister(rDCM.FIXEDSEED))

    Ep_A = [-0.5231677865060994 0.0 -0.7986242409570947;
     -1.1850222856258559 -0.49886852758295025 0.0;
      0.0 0.0 -0.4259224542476115]
    Ep_C = [-1.14490153172882 0.15614346264074028 -2.641991008076796 0.18702790710363;
      0.0 0.0 1.0033099014594844 0.5181487878771377;
        -0.46860588216767457 0.0 1.0823812056084292 0.0]

    Ep_transit = zeros(3)
    Ep_decay   = zeros(3)
    Ep_epsilon = 0.0

    @test all(Ep_A .≈ A)
    @test all(Ep_C .≈ C)
    @test all(Ep_transit .≈ transit)
    @test all(Ep_decay .≈ decay)
    @test Ep_epsilon ≈ epsilon
  end
end

function test_get_rigid_prior(dcm)

  @testset "Rigid prior" begin
    rdcm = RigidRdcm(dcm)
    μ0, l0, a0, β0 = rDCM.get_priors(rdcm)

    nr = size(dcm.a,1)

    pE_A = [-0.500000000000000	0	0	0	0	0	0	0	0	0;
    0	-0.500000000000000	0	0	0	0	0	0	0	0;
    0	0	-0.500000000000000	0	0	0	0	0	0	0;
    0	0	0	-0.500000000000000	0	0	0	0	0	0;
    0	0	0	0	-0.500000000000000	0	0	0	0	0;
    0	0	0	0	0	-0.500000000000000	0	0	0	0]

    pE_C = [0	0	0	0	0	0;
    0	0	0	0	0	0;
    0	0	0	0	0	0;
    0	0	0	0	0	0;
    0	0	0	0	0	0;
    0	0	0	0	0	0;
    0	0	0	0	0	0]

    pC_A = [400	6.25000000000000	Inf	Inf	6.25000000000000	Inf	Inf	Inf	Inf	Inf;
    Inf	400	6.25000000000000	Inf	Inf	Inf	Inf	Inf	Inf	Inf;
    Inf	Inf	400	6.25000000000000	Inf	Inf	Inf	6.25000000000000	Inf	Inf;
    Inf	Inf	Inf	400	6.25000000000000	Inf	Inf	Inf	Inf	Inf;
    Inf	Inf	Inf	Inf	400	Inf	Inf	Inf	Inf	Inf;
    Inf	Inf	Inf	Inf	Inf	400	6.25000000000000	Inf	Inf	6.25000000000000]

    pC_C = [Inf	1	Inf	Inf	Inf	Inf;
    Inf	Inf	Inf	Inf	Inf	Inf;
    Inf	Inf	1	Inf	Inf	Inf;
    Inf	Inf	Inf	Inf	Inf	Inf;
    Inf	Inf	Inf	Inf	1	Inf;
    Inf	Inf	Inf	Inf	Inf	Inf;
    Inf	Inf	Inf	Inf	Inf	1]

    @test all(pE_A .≈ μ0[1:6,1:10])
    @test all(pE_C .≈ μ0[1:7,nr+1:nr+6])
    @test all(pC_A .≈ l0[1:6,1:10])
    @test all(pC_C .≈ l0[1:7,nr+1:nr+6])

    @test a0 == 2.0
    @test β0 == 1.0
  end
end

function test_prior()

  dcm = load_example_DCM()

  a = BitArray([1 0 1;1 1 0;0 0 1])
  c = BitArray([1 1 1 1;0 0 1 1;1 0 1 0])

  @testset verbose=true "Priors" begin
    test_pE(a,c)
    test_sample_prior(a,c)
    test_get_rigid_prior(dcm)
  end

end

test_prior()
