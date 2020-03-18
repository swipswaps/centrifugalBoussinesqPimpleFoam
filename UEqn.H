    // Solve the momentum equation

	vector centrifugal(0,1,0); 

	
	// Now adjust the nudging. 
	forAll(NudgingTerm.internalField(), celli) {
		label centercelli =  CenterLookup[celli]; 
		
		if (centercelli == -1) continue; // ignore phantom cells because of the cyclic/domain decompostion boundary cells. 
		NudgingTerm.internalField()[celli].component(0) = nudgingFactor*(U.internalField()[centercelli].component(0)- Umean.internalField()[centercelli].component(0));
		
	}
	
	forAll(NudgingTerm.boundaryField(), celli)
	{
		NudgingTerm[celli].component(0) = 0; 
		NudgingTerm[celli].component(1) = 0; 
		NudgingTerm[celli].component(2) = 0; 
	}	
	
	
	
    fvVectorMatrix UEqn(
		(nonlinear ? fvm::ddt(U) + fvm::div(phi, U) : fvm::ddt(U))
	      + turbulence->divDevReff(U)
//	      - CentrifugalFactor*centrifugal*sqr(U.component(0))/Rcurvature	      
	      + NudgingTerm
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


