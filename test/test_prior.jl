
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
    A, C, transit, decay, epsilon = sample_from_prior(a,c,prior;rng=Xoshiro(rDCM.FIXEDSEED))

    Ep_A = [-0.4671518499331545 0.0 -1.9553465718263805;
     -2.3462895892117315 -0.5364080617587345 0.0;
      0.0 0.0 -0.52925813328157]
    Ep_C = [-0.07258186804586991 0.6316208311167526 1.4386832757114134 -1.175371133217587;
     0.0 0.0 -0.35416984337044605 -0.5933950393067663;
      -0.9129233863399265 0.0 0.796126919278033 0.0]


    Ep_transit = zeros(3)
    Ep_decay   = zeros(3)
    Ep_epsilon = 0.0

    @test all(Ep_A .≈ A)
    @test all(Ep_C .≈ C)
    @test all(Ep_transit .≈ transit)
    @test all(Ep_decay .≈ decay)
    @test Ep_epsilon ≈ epsilon

    # sample also HRF params
    transit_ref = [0.021388841464128085, 0.0474136423020816, 0.006429462759905607]
    decay_ref = [-0.01467592799178971, -0.018633728172154647, 0.058226560968698854]
    epsilon_ref = 0.0006366637688044764
    _, _, transit, decay, epsilon = sample_from_prior(a, c; fixHRF=false, rng=Xoshiro(rDCM.FIXEDSEED))

    @test all(transit_ref .≈ transit)
    @test all(decay_ref .≈ decay)
    @test epsilon_ref ≈ epsilon
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
