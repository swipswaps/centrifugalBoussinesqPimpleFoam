    singlePhaseTransportModel laminarTransport(U, phi);

    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta(laminarTransport.lookup("beta"));

    //Info << "Reading Anisotropic diffusion " << endl;
    //dimensionedTensor AnistropicDiffusion(laminarTransport.lookup("anisotropicdiffusion"));
    //Info << "Reading SpongeLayer diffusion " << endl;
    //dimensionedTensor SpongeAnistropicDiffusion(laminarTransport.lookupOrDefault("Spongeanisotropicdiffusion",AnistropicDiffusion));
    //dimensionedScalar SpongeLayerThickness = laminarTransport.lookupOrDefault("SpongeThickness",dimensionedScalar("st",dimLength,scalar(0)));
    //scalar            TransitionSpongeLayerFraction = laminarTransport.lookupOrDefault("TransitionSpongeLayerFraction",0.5);

    // Reference temperature [K]
	// dimensionedScalar TRef(laminarTransport.lookup("TRef"));

    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

	// read radius of curvature. 
	dimensionedScalar Rcurvature(laminarTransport.lookup("Rcurvature"));

	
