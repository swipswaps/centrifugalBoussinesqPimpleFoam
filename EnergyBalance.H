




/*{
		Info << "Momentum balance - all domain" << endl;
		Info << "-----------------------------" << endl;
                Info << fvc::ddt(U)->weightedAverage(mesh.V()) << endl; 
                Info << fvc::div(phi, U)->weightedAverage(mesh.V()) << endl; 
                Info << fvc::laplacian(AnisotropicDiffusion,U)->weightedAverage(mesh.V()) << endl;
                Info << fvc::grad(p_rgh)->weightedAverage(mesh.V()) << endl;
                Info << (g*rhok_tag)->weightedAverage(mesh.V()) << endl;
		Info << NudgingTerm.weightedAverage(mesh.V()) << endl;

		Info << (
                (nonlinear ? fvc::ddt(U) + fvc::div(phi, U) : fvc::ddt(U))
              + NudgingTerm
              - fvc::laplacian(AnisotropicDiffusion,U)
              + fvc::grad(p_rgh)
              - g*rhok_tag)->weightedAverage(mesh.V()) << endl;



		Info << "Momentum balance - ZONE " << endl;
		Info << "-----------------------------" << endl;
                Info << "ddt " << IntegrateVector( fvc::ddt(U)*zoneSelector)                       		 / ZoneVolume  << endl; 
                Info << "div(phi,U) " << IntegrateVector( fvc::div(phi, U)*zoneSelector)                       / ZoneVolume  << endl; 
                Info << "laplcian "   << IntegrateVector( fvc::laplacian(AnisotropicDiffusion,U)*zoneSelector) / ZoneVolume  << endl; 
                Info << "grad "       << IntegrateVector( fvc::grad(p_rgh)*zoneSelector) / ZoneVolume  << endl; 
                Info << "grhok "      << IntegrateVector( (g*rhok_tag)*zoneSelector) / ZoneVolume  << endl; 
		Info << "Nudging "    << IntegrateVector( NudgingTerm*zoneSelector) / ZoneVolume  << endl;  

		Info << " All " <<  IntegrateVector( (
                (nonlinear ? fvc::ddt(U) + fvc::div(phi, U) : fvc::ddt(U))
              + NudgingTerm
              - fvc::laplacian(AnisotropicDiffusion,U)
              + fvc::grad(p_rgh)
              - g*rhok_tag)*zoneSelector) / ZoneVolume  << endl;  


		
		Info << "Energy balance - all domain" << endl;
		Info << "-----------------------------" << endl;
                Info << (U&fvc::ddt(U))->weightedAverage(mesh.V()) << endl; 
                Info << (U&fvc::div(phi, U))->weightedAverage(mesh.V()) << " || " << (fvc::div(phi,0.5*U&U))->weightedAverage(mesh.V())  <<endl; 
                Info << (U&fvc::laplacian(AnisotropicDiffusion,U))->weightedAverage(mesh.V()) << endl;
                Info << fvc::div(phi,p_rgh)->weightedAverage(mesh.V()) << endl;
                Info << (U&(g*rhok_tag))->weightedAverage(mesh.V()) << endl;
		Info << (U&NudgingTerm)->weightedAverage(mesh.V()) << endl;

		Info << (U&(
                (nonlinear ? fvc::ddt(U) + fvc::div(phi, U) : fvc::ddt(U))
              + NudgingTerm
              - fvc::laplacian(AnisotropicDiffusion,U)
              + fvc::grad(p_rgh)
              - g*rhok_tag))->weightedAverage(mesh.V()) << endl;

}*/
