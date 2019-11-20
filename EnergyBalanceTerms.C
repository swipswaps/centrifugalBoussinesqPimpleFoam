#include "EnergyBalanceTerms.H" 

EnergyBalanceTerms::EnergyBalanceTerms(
				fvMesh& mesh,
				Time&   runTime,
				volVectorField &U,
				surfaceScalarField &phi,
				volTensorField &AnisotropicDiffusion,
				volScalarField &p_rgh,
				volScalarField &T,
				word ZoneName) :
	_calcMean(false),
	_work(true),
	mesh(mesh),
	runTime(runTime),
	U(U), 
	phi(phi), 
	AnisotropicDiffusion(AnisotropicDiffusion), 
	p_rgh(p_rgh), 
	T(T), 
	zoneSelector(
            IOobject
            (
                    "ZoneSelector",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
	),
	mean_energy_dUdt(0)
{

	IOobject check_mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE);
	if (check_mean_U_Header.headerOk()) { 
		_calcMean = false;
		IOobject mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		IOobject mean_p_rgh_Header("mean_p_rgh" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		IOobject mean_T_Header("mean_T" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		IOobject mean_Phi_header("mean_phi",mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);

		mean_U 		= new Foam::volVectorField(mean_U_Header,mesh); 
		mean_p_rgh  	= new Foam::volScalarField(mean_p_rgh_Header,mesh,dimensionedScalar("uu",p_rgh.dimensions(),0)); ;
		mean_T  	= new Foam::volScalarField(mean_U_Header,mesh,dimensionedScalar("uu",T.dimensions(),0)); ;
		mean_phi	= new Foam::surfaceScalarField(mean_U_Header,mesh,dimensionedScalar("uu",phi.dimensions(),0)); 

	} else {
		Info << "calculate mean " << endl;
		if (ZoneName != "") { 
			_work = false; 
		} else {
			_calcMean = true;

			IOobject mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);
			IOobject mean_p_rgh_Header("mean_p_rgh" ,mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);
			IOobject mean_T_Header("mean_T" ,mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);
			IOobject mean_Phi_Header("mean_phi",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);

			mean_U 		= new Foam::volVectorField(mean_U_Header,mesh,dimensionedVector("uu",U.dimensions(),vector::zero)); 
			mean_p_rgh  	= new Foam::volScalarField(mean_p_rgh_Header,mesh,dimensionedScalar("uu",p_rgh.dimensions(),0)); ;
			mean_T  	= new Foam::volScalarField(mean_T_Header,mesh,dimensionedScalar("uu",T.dimensions(),0)); ;
			mean_phi	= new Foam::surfaceScalarField(mean_Phi_Header,mesh,dimensionedScalar("uu",phi.dimensions(),0)); 
		}
	}

	// Specialized code for the nudging. 
    	meshSearch meshSearch(mesh); 
	scalar CenterOfDomain = (max(mesh.C().component(1))-min(mesh.C().component(1))).value()/2;    

	
/*
	label cellzoneID  =  mesh.cellZones().findZoneID("Work");
	const labelList& cellzonelist =  mesh.cellZones()[cellzoneID];
	scalar tsize = cellzonelist.size() ;
	reduce (tsize,sumOp<scalar>()); 

	Info << "zone work has " << tsize << " cells " << endl;
	forAll(cellzonelist,cellindx) {
		label currentZoneIndx = cellzonelist[cellindx];
	
		if (globalmeshdata.parallel()) {
		
			label lbl = -1;

			if (globalindexer.isLocal(currentZoneIndx)) {
				 lbl = globalindexer.toLocal(currentZoneIndx);
	 			 zoneSelector[lbl] = mesh.V()[currentZoneIndx];	
			}
		}
		else { 
			zoneSelector[currentZoneIndx] = mesh.V()[currentZoneIndx];
		}
	}
*/

}


void EnergyBalanceTerms::update() {
	if (_work) { 
		if (_calcMean) { 
			update_mean();
		}
	}
}


void EnergyBalanceTerms::finalize() { 
	if (_work) { 
		if (_calcMean) { 
			finalize_calculate_mean();

			// now should write in the first time step. 
//			const instantList& times  = mesh.time().times();


			runTime.setTime(mesh.time().startTime(),0); 

			mean_U->write();
			mean_p_rgh->write();
			mean_phi->write();
			mean_T->write(); 

		}
	}
}


vector EnergyBalanceTerms::momentum_dUdt() {
	return Integrate(fvc::ddt(U));
}

scalar EnergyBalanceTerms::energy_dUdt() { 
	scalar dt = mesh.time().deltaTValue();

	scalar local_balance = Integrate( U&fvc::ddt(U) )*dt;
	mean_energy_dUdt += local_balance;

	return local_balance;
}

	// ====================================== Convection ====================	
//- Calculates the <div(phi,U)>
vector EnergyBalanceTerms::momentum_Convection() {
	return Integrate(fvc::div(phi,U)); 
}

//- Calculates the <U&div(phi,U)> and accumulates it for the average. 
scalar EnergyBalanceTerms::energy_Convection() {
	
	scalar dt = mesh.time().deltaTValue();

	scalar local_balance = Integrate( U&fvc::div(phi,U) )*dt;
	mean_energy_Convection += local_balance;

	return local_balance;
}

	// ====================================== Diffusion =====================
//- Calculates the <laplacian(AnisitropicDiffusion,U)>
vector EnergyBalanceTerms::momentum_Diffusion() {
	return Integrate(fvc::laplacian(AnisotropicDiffusion,U)); 
}

//- Calculates the <U&laplacian(AnisitropicDiffusion,U)> and accumulates it for the average. 
scalar EnergyBalanceTerms::energy_Diffusion() {

	scalar dt = mesh.time().deltaTValue();

	scalar local_balance = Integrate( U&fvc::laplacian(AnisotropicDiffusion,U) )*dt;
	mean_energy_Diffusion += local_balance;

	return local_balance;
}

	// ====================================== Pressure  =====================
//- Calculates the <grad(p_rgh)>
vector EnergyBalanceTerms::momentum_Pressure() {
	return Integrate(fvc::grad(p_rgh)); 
}

//- Calculates the <U&grad(p_rgh)> and accumulates it for the average. 
scalar EnergyBalanceTerms::energy_Pressure() {
	scalar dt = mesh.time().deltaTValue();

	scalar local_balance = Integrate( U&fvc::grad(p_rgh) )*dt;
	mean_energy_Diffusion += local_balance;

	return local_balance;
} 

	// ====================================== Nudging =======================





	// ==========================================================
	// ==========================================================
void EnergyBalanceTerms::update_mean() { 

	scalar dt = mesh.time().deltaTValue();


	volVectorField& meanU 		= *mean_U;
	volScalarField& meanP 		= *mean_p_rgh;
	surfaceScalarField& meanphi 	= *mean_phi; 
	volScalarField& meanT  		= *mean_T; 

	meanU   += U*dt; 
	meanP   += p_rgh*dt; 
	meanphi += phi*dt; 
	meanT   += T*dt;
	
	
}


void EnergyBalanceTerms::finalize_calculate_mean() {
	scalar TotalTime = (mesh.time().endTime()-mesh.time().startTime()).value();

	Info << " Total running time is " << TotalTime << endl;
	volVectorField& meanU 		= *mean_U;
	volScalarField& meanP 		= *mean_p_rgh;
	surfaceScalarField& meanphi 	= *mean_phi; 
	volScalarField& meanT  		= *mean_T; 

	meanU   /= TotalTime; 
	meanP   /= TotalTime; 
	meanphi /= TotalTime; 
	meanT   /= TotalTime;
	
}


