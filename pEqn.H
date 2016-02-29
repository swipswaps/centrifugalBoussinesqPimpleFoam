{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    volScalarField g_gradptag = (g&fvc::grad(rAU*rhok_tag));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + rAUf*fvc::ddtCorr(U, phi)
      );


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA) + g_gradptag 
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
	    /**

		Note that it should be g&mesh.Sf() because we calculated d/dz.  = g*Sz
		It is NOT equivalent to g&(0,0,1)*magSf(). because it is sqrt(Sx**2+Sy**2+Sz**2). 
		and therefore, g&(0,0,1)*magSf() = 9.8*sqrt(Sx**2+Sy**2+Sz**2). 
	    */
            phi = phiHbyA - p_rghEqn.flux() + (g&mesh.Sf())*fvc::interpolate(rAU*rhok_tag); // remember that g is negative!

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((-p_rghEqn.flux())/rAUf) + rAU*g*rhok_tag; // remember that g is negative!
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"
    p = p_rgh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p; 
    }
}
