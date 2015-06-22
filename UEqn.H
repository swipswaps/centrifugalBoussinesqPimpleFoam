    // Solve the momentum equation

	vector centrifugal(0,1,0); 

	
	
    fvVectorMatrix UEqn(
		(nonlinear ? fvm::ddt(U) + fvm::div(phi, U) : fvm::ddt(U))
	      + turbulence->divDevReff(U)
	      - CentrifugalFactor*centrifugal*sqr(U.component(0))/Rcurvature	      
	      + nudgingFactor*(U-Umean)
	      - fvm::laplacian(AnisotropicDiffusion,U)
	     ==
		fvOptions(U)
    );


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


