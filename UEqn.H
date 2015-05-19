    // Solve the momentum equation

	vector centrifugal(0,1,0); 

	
	
	// Calculating the nudging term
	forAll(NudgingTerm.internalField(), celli) {
		
		vector position = mesh.C().internalField()[celli];
		//Info << "For the cell in position " << position;
		position.component(1) = CenterOfDomain;
		//Info << " find the cell in  " << position;
		label centercelli = mesh.findCell(position);
		//Info << " --> " << centercelli << " || (U-Umean)=" << U.component(0)->operator[](centercelli) << " - " <<  Umean.component(0)->operator[](centercelli) << endl;
		NudgingTerm.internalField()[celli].component(0) = nudgingFactor*(U.component(0)->operator[](centercelli) - Umean.component(0)->operator[](centercelli));
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
	      - CentrifugalFactor*centrifugal*sqr(U.component(0))/Rcurvature	      
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


