    // Solve the momentum equation

	vector centrifugal(0,1,0); 

    fvVectorMatrix UEqn(
		(nonlinear ? fvm::ddt(U) + fvm::div(phi, U) : fvm::ddt(U))
	      + fvm::div(phimean, U)
	      + fvc::div(phi, Umean)
	      + turbulence->divDevReff(U)
	      - CentrifugalFactor*centrifugal*sqr(Umean.component(0)+U.component(0))/Rcurvature	      
	      - fvm::laplacian(AnisotropicDiffusion,U)
	     ==
		fvOptions(U)
    );

//+U.component(0)

    //if (whiteNoiseFlag && !ExplicitwhiteNoiseFlag) { 
//	  UEqn.source() += Uwhitenoise;
    //} 

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
	  - fvc::grad(p_rgh)
          + g*rhok_tag
        );

        fvOptions.correct(U);
    }


