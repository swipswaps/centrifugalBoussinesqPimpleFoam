/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
	centrifugalBoussinesqPimpleFoam

Description
	* Takes into account the centrifugal term. 
	* The Umean is a nudging term. we use the nudging term as with the center of the cell. 
	* 
 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
//#include "RASModel.H"
#include "radiationModel.H"
#include "turbulentTransportModel.H"

#include "fvIOoptionList.H"
#include "pimpleControl.H"
#include "interpolation.H"
#include "Random.H"
#include "meshSearch.H"

#include <sstream>
#include <fstream>

#include "globalMeshData.H"
#include "globalIndex.H"

#include "EnergyBalanceTerms.H" 


//#include "stdlib.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "EnergyFunctions.H"
int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
//    #include "createIncompressibleRadiationModel.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    Info << " \n\nCentrifugal solver  :  0.2.3 - with balances" << endl; 
    Info << " -------------------------- " << endl; 

    const bool nonlinear = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("nonlinearSolver", true);

    const bool    ExplicitwhiteNoiseFlag       = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("explicitwhitenoise", true);
    const bool    whiteNoiseFlag       	       = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("whitenoise", false);
    const scalar  whiteNoiseSeed               = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("whitenoise_seed", 0);
    
//    const dimensionedScalar nudgingFactor(mesh.solutionDict().subDict("PIMPLE").lookup("nudging"));
    const scalar nudgingFactor				   = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("nudging",0.0);

	
    IOdictionary MeanKineticEnergy
    (
	IOobject
	(
		"MeanKineticEnergy",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    	    )
     );
     dimensionedScalar KE(MeanKineticEnergy.lookup("MeanKineticEnergy"));
     scalar EnergyFraction(MeanKineticEnergy.lookupOrDefault("EnergyFraction",0.01));

     dimensionedScalar rootdT = sqrt(runTime.deltaT());

     scalar KEfactor(EnergyFraction*KE.value()/rootdT.value());
    Info << "* " << ( nonlinear ? "non-linear" : "linear") << " solver. --- " <<endl; 

    Info << "\n ---- Running " << ( nonlinear ? "non-linear" : "linear") << " solver. --- " <<endl; 
    if (whiteNoiseFlag) { 
	Info << "\n ---- Using whitenoise, seed : " << whiteNoiseSeed << endl; 
	Info << "\t ---- KineticEnertgy  " << KE.value()  << " | Fraction " << EnergyFraction << endl; 
    } else { 
	    Info << "* no white noise " << endl; 
    }
    Info << "Nudging coefficient " << nudgingFactor << endl; 

    label n=Pstream::nProcs();
    Random perturbation(whiteNoiseSeed+Pstream::myProcNo()*n);	



//    #include "EnergyInit.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    EnergyBalanceTerms domainEnergyBalanceTerms(mesh,runTime,U,phi,AnisotropicDiffusion,p_rgh,T,"");



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
	scalar MinDiffusionX = min(AnisotropicDiffusion.component(0)).value(); 
	scalar MaxDiffusionX = max(AnisotropicDiffusion.component(0)).value(); 		
	
	scalar MinDiffusionY = min(AnisotropicDiffusion.component(4)).value(); 
	scalar MaxDiffusionY = max(AnisotropicDiffusion.component(4)).value(); 		

	scalar MinDiffusionZ = min(AnisotropicDiffusion.component(8)).value(); 
	scalar MaxDiffusionZ = max(AnisotropicDiffusion.component(8)).value(); 		
	
	
	const scalar CenterOfDomain = (max(mesh.C().component(1))-min(mesh.C().component(1))).value()/2;
	
	List<int> CenterLookup(mesh.C().internalField().size());

	Info << " Building the lookup table for the center  " << CenterOfDomain  << " mesh bounds [" << max(mesh.C().component(1)).value() << " , " << min(mesh.C().component(1)).value() << "]" << endl; 
 
	meshSearch searchEngine(mesh);
	label centercelli = -1;
	forAll(mesh.C().internalField(), celli) { 
		if (celli % 10000 == 0) {
			Info << celli << "/" << mesh.C().internalField().size() << endl;
		}
		vector position = mesh.C().internalField()[celli];
		position.component(1) = CenterOfDomain;
		centercelli = searchEngine.findCell(position,centercelli); //findCell(position);
		
		if (centercelli == -1) { 
			Sout << " Celll not found  in processor "  << Pstream::myProcNo() << ": At cell " 
			     << mesh.C().internalField()[celli] << " label " << celli << " looking for position " << position << " Not Found!" << endl;
		}
		CenterLookup[celli] = centercelli; 
		
	}
	Info << " -- End -- " << endl; 

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

		if (whiteNoiseFlag) {
		
		
		  forAll(mesh.C().internalField(), celli)
		  {
		  
			scalar factorX; 
			scalar factorY;
			scalar factorZ;
			scalar factorT;
			
			if (MaxDiffusionX-MinDiffusionX < 1e-10) { 
				factorX = 1;
				factorT = 1;
				
			} else {
				factorX = ((MaxDiffusionX-AnisotropicDiffusion[celli].component(0))/(MaxDiffusionX-MinDiffusionX) );
				factorT = ((MaxDiffusionX-AnisotropicDiffusion[celli].component(0))/(MaxDiffusionX-MinDiffusionX) );
			}

			if (MaxDiffusionY-MinDiffusionY < 1e-10) { 
				factorY = 1;
			} else {
				factorY = ((MaxDiffusionY-AnisotropicDiffusion[celli].component(4))/(MaxDiffusionY-MinDiffusionY) );
			}
			
			if (MaxDiffusionZ-MinDiffusionZ < 1e-10) { 
				factorZ = 1;
			} else {
				factorZ = ((MaxDiffusionZ-AnisotropicDiffusion[celli].component(8))/(MaxDiffusionZ-MinDiffusionZ) );
			}
		  
			Uwhitenoise[celli].component(0) = KEfactor*perturbation.GaussNormal()*factorX; 
			Uwhitenoise[celli].component(1) = KEfactor*perturbation.GaussNormal()*factorY; 
			Uwhitenoise[celli].component(2) = KEfactor*perturbation.GaussNormal()*factorZ; 
			Twhitenoise[celli]              = KEfactor*perturbation.GaussNormal()*factorT; 
			
			
		  }
		  
		  forAll(mesh.C().boundaryField(), celli)
		  {
			Uwhitenoise[celli].component(0) = 0; 
			Uwhitenoise[celli].component(1) = 0; 
			Uwhitenoise[celli].component(2) = 0; 
			Twhitenoise[celli]              = 0; 
		  }
		
		}
	


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
		
          #include "UEqn.H"
	  #include "TEqn.H"
	    
	        // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
	
		if (whiteNoiseFlag && ExplicitwhiteNoiseFlag) { 
			U +=  Uwhitenoise*runTime.deltaT()/dimensionedScalar("corrector",dimTime,scalar(1));
			T +=  Twhitenoise*runTime.deltaT()/dimensionedScalar("corrector",dimTime,scalar(1));
		}

		
		domainEnergyBalanceTerms.update();
		//#include "EnergyBalance.H"

	{

	wordList direction(3);
	direction[0] = "x";
	direction[1] = "y";
	direction[2] = "z";


	wordList component(3);
	component[0] = "u";
	component[1] = "v";
	component[2] = "w";

	
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


		Info << "Total balance " <<  (fvc::ddt(U) + fvc::div(phi, U)+ NudgingTerm- fvc::laplacian(AnisotropicDiffusion,U)  + fvc::grad(p_rgh) - g*rhok_tag)->weightedAverage(mesh.V()) << endl;
		Info << "Total Enegy balance " <<  (U&(fvc::ddt(U) + fvc::div(phi, U)+ NudgingTerm- fvc::laplacian(AnisotropicDiffusion,U)  + fvc::grad(p_rgh) - g*rhok_tag))->weightedAverage(mesh.V()) << endl;

		dimensionedScalar dt_nonlin = (U&(fvc::ddt(U) + fvc::div(phi, U)))->weightedAverage(mesh.V());
		Info << "--- Term wise --- " << endl;
		Info << "dt  " << (U&fvc::ddt(U))->weightedAverage(mesh.V()) <<endl;
		Info << "convection  " << (U&fvc::div(phi, U))->weightedAverage(mesh.V()) <<endl;
		Info << "\t convection2 " << (U&fvc::div(U*U))->weightedAverage(mesh.V()) << endl;
		Info << "\t convection3 " << (fvc::div(phi,0.5*U& U))->weightedAverage(mesh.V()) << endl;
		Info << "\t convection4 " << (fvc::div(U*(0.5*U& U)))->weightedAverage(mesh.V()) << endl;
		Info << "Nudging  " << (U&NudgingTerm)->weightedAverage(mesh.V()) <<endl;
		Info << "Laplacian " << -(U&fvc::laplacian(AnisotropicDiffusion,U))->weightedAverage(mesh.V())<< endl;
		Info << "pressure " << (U&fvc::grad(p_rgh))->weightedAverage(mesh.V()) << endl;
		Info << "buoyancy " << (-U&g*rhok_tag)->weightedAverage(mesh.V())<< endl;

		Info << "--- Cell wise --- " << endl; 
		volVectorField phiU(IOobject("phiU",mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),fvc::div(phi,U)); 
		volVectorField divUU(IOobject("divUU",mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),fvc::div(phi*fvc::interpolate(U)));
		volTensorField divUUterms(IOobject("gradUU",mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),fvc::grad(phi*fvc::interpolate(U)));  
		volTensorField uu = U*U;
		phiU.write();
		divUU.write();
		
		Info << "--- Momentum --- " << endl; 
		Info << (fvc::div(phi, U))->weightedAverage(mesh.V()) <<endl;
		Info << (fvc::div(U*U))->weightedAverage(mesh.V()) <<endl;

		Info << (fvc::div(phi, U))->component(0)->weightedAverage(mesh.V()) <<endl;
		Info << (fvc::div(U*U)   )->component(0)->weightedAverage(mesh.V()) <<endl;

		volTensorField gradUU = fvc::grad(U);


		int global=0;
		for (int i=0 ; i<3 ;i++) {
			for (int j=0 ; j<3 ; j++) { 
				gradUU.component(global) = fvc::grad(uu.component(global))->component(j);
				global++; 
			}
		}


		forAll(mesh.C(),celli) { 
			Info << mesh.C()[celli] << ",";
			Info << mesh.V()[celli] << ",";
			Info << celli << "," << phiU[celli].component(0) << "," <<divUU[celli].component(0) << "," ;
			for (int i=0 ; i < 9 ; i++ ) {
				Info << divUUterms[celli].component(i); 
				if (i<8) Info << ","; 
			}
			Info << endl;


		//	Info << " div(phi,U)->component(0) = d(uu)/dx+d(vu)/dy+d(wu)/dz =  " << phiU[celli].component(0) << " || ";
		//	Info << "fvc::grad(U*U): d(uu)/dx+d(vu)/dy+d(wu)/dz =" << gradUU[celli].component(0) << "+" <<  gradUU[celli].component(1) << "+" << gradUU[celli].component(2) << " = " << 
		//			gradUU[celli].component(0) + gradUU[celli].component(1) + gradUU[celli].component(2) << endl;
		}

		label cellid = 23400;
		// Now calculate the divergance for cell 23400 by hand using the gauss theorem. 

		// Then calculate the grad for cell 23400  by hand using the gauss theorem. 


		



		/*
		volVectorField Ekcomponents = U;
		Ekcomponents.component(0) = 0.5*uu.component(0);
		Ekcomponents.component(1) = 0.5*uu.component(4);
		Ekcomponents.component(2) = 0.5*uu.component(8);


		Info << "--- energy --- " << endl; 
		for (int i=0 ; i<3 ;i++) {
			Info << "All convection  " << (U.component(i)*fvc::div(phi, U)->component(i))->weightedAverage(mesh.V()) <<endl;
			Info << "\t  " << (U.component(i)*fvc::div(phi*fvc::interpolate(U))->component(i))->weightedAverage(mesh.V()) <<endl;
			for (int j=0 ; j<3 ; j++) { 
				Info << "\t\t "<< direction[i] << direction[j] << (U.component(i)*fvc::grad(phi*fvc::interpolate(U))->component(j))->weightedAverage(mesh.V()) <<endl;
			}

		}


		Info << "Ek " << (fvc::div(U*Ekcomponents))->weightedAverage(mesh.V()) <<endl;

		int global=0;
		for (int i=0 ; i<3 ;i++) {
			for (int j=0 ; j<3 ; j++) { 
				Info << "Momentum " << direction[i] << direction[j] << ": " << fvc::grad(uu.component(global))->component(i)->weightedAverage(mesh.V()) << endl;
				global++; 

			}


		}
		Info << "\t\t -- " <<endl;
		global=0;
		for (int i=0 ; i<3 ;i++) 
			for (int j=0 ; j<3 ; j++) { 
				Info << "Energy " << direction[i] << direction[j] << ": " << (U.component(j)*fvc::grad(uu.component(global))->component(i))->weightedAverage(mesh.V()) << " || " << endl;
				Info << "\t" << fvc::grad((U*Ekcomponents)->component(global))->component(i)->weightedAverage(mesh.V()) << endl;
				global++; 

		}
				

		Info << (U&(fvc::div(phi, U)))->weightedAverage(mesh.V()) <<endl;
		Info << (U&(fvc::div(U*U)))->weightedAverage(mesh.V()) <<endl;
*/
	}



		runTime.write();
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
            
    }


    domainEnergyBalanceTerms.finalize();
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
//		outname << "rand" << Pstream::myProcNo();
//		std::ofstream out;
		//out.open(outname.str().c_str()); 
		
		//for (long jj=0;jj<10000;jj++){ 
//			out << perturbation.GaussNormal() << std::endl;
		//}

		
		// Get the random field. 
		//if (Pstream::myProcNo() == 0)
		//{
			//List<scalar> allV(n);
			//for(label i=1; i<n; i++)
			//{
				// create the input stream from processor i
				//IPstream vStream(Pstream::blocking, i);
				//vStream >> allV[i];
			//}
			//Info << allV << endl;
			//exit(1);
		//} else { 
			
						
			//OPstream vectorStream(Pstream::blocking, 0);
			//vectorStream << perturbation.GaussNormal();
			
//		}
