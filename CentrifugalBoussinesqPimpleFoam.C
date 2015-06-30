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
#include "RASModel.H"
#include "radiationModel.H"
#include "fvIOoptionList.H"
#include "pimpleControl.H"
#include "interpolation.H"
#include "Random.H"
#include "meshSearch.H"

#include <sstream>
#include <fstream>



//#include "stdlib.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createIncompressibleRadiationModel.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    Info << " \n\nCentrifugal solver  :  0.2.3" << endl; 
    Info << " -------------------------- " << endl; 
    Info << " \n\nBased on basestate version 0.8.0" << endl; 

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



    pimpleControl pimple(mesh);

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
			Sout << " Celll not found  in processor "  << Pstream::myProcNo() << ": At cell " << mesh.C().internalField()[celli] << " label " << celli << " looking for position " << position << " Not Found!" << endl;
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

		Utotal = U; 
		Ttotal = T+Tmean;

		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
            
    }

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
