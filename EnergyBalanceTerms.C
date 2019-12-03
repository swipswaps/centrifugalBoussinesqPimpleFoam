#include "EnergyBalanceTerms.H" 


#define zeroScalar(variable,dims) \
	variable(		\
            IOobject		\
            (			\
                    "##variable##",		\
                    mesh.time().timeName(),	\
                    mesh,			\
                    IOobject::NO_READ,		\
                    IOobject::AUTO_WRITE	\
            ),					\
            mesh,				\
            dimensionedScalar("zero",dims,0.0)\
	)

#define oneScalar(variable,dims) \
	variable(		\
            IOobject		\
            (			\
                    "##variable##",		\
                    mesh.time().timeName(),	\
                    mesh,			\
                    IOobject::NO_READ,		\
                    IOobject::AUTO_WRITE	\
            ),					\
            mesh,				\
            dimensionedScalar("zero",dims,1.0)\
	)

#define zeroVector(variable,dims) \
	variable(		\
            IOobject		\
            (			\
                    "##variable##",		\
                    mesh.time().timeName(),	\
                    mesh,			\
                    IOobject::NO_READ,		\
                    IOobject::AUTO_WRITE	\
            ),					\
            mesh,				\
            dimensionedVector("zero",dims,vector::zero)\
	)

#define zeroTensor(variable,dims) \
	variable(		\
            IOobject		\
            (			\
                    "##variable##",		\
                    mesh.time().timeName(),	\
                    mesh,			\
                    IOobject::NO_READ,		\
                    IOobject::AUTO_WRITE	\
            ),					\
            mesh,				\
            dimensionedTensor("zero",dims,tensor::zero)\
	)

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
	zeroScalar(zoneSelector,dimless),
	zeroVector(dmeanUdt,U.dimensions()),
	zeroVector(tagU,U.dimensions()),
	zeroVector(tagU_ls,U.dimensions()),
	zeroScalar(tag_p_rgh,p_rgh.dimensions()),
	zeroScalar(tag_T,T.dimensions()),
	zeroTensor(reynoldsU,U.dimensions()*U.dimensions())
{

	energy_variables_set_zero(dUdt); 
	energy_variables_set_zero(pressure);


	IOobject check_mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE);
	if (check_mean_U_Header.headerOk()) { 
		Info << "=========================== Mean exists: loading " << endl;
		_calcMean = false;
		IOobject mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		IOobject mean_p_rgh_Header("mean_p_rgh" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		IOobject mean_T_Header("mean_T" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		IOobject mean_Phi_header("mean_phi",mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);

		mean_U 		= new Foam::volVectorField(mean_U_Header,mesh); 
		mean_p_rgh  	= new Foam::volScalarField(mean_p_rgh_Header,mesh);
		mean_T  	= new Foam::volScalarField(mean_T_Header,mesh);
		mean_phi	= new Foam::surfaceScalarField(mean_Phi_header,mesh); 

		Info << " -- calculating the first u' " << endl;
		tagU = U-(*mean_U); 

		Info << " Loading the last time step to calculate d[ubar]/dt " << endl;
		word lastTime(name(runTime.endTime().value()));
		dimensionedScalar timeSpan = runTime.endTime() - runTime.startTime();

		// get the dubar/dt = (u|1-u|0)/Time
		Info << " Reading time step " << lastTime << endl;
		volVectorField Ulast(IOobject("U",lastTime,mesh,IOobject::MUST_READ,IOobject::NO_WRITE),
				      mesh);


		// calculate the mean fields. 
		volVectorField& meanU 		= *mean_U;

		mean_energy_dUdt 	= Integrate( (meanU) & ((Ulast-U)/timeSpan)   );
		mean_energy_pressure 	= Integrate( (meanU) & fvc::grad(*mean_p_rgh) );


		// Tests for the later calculation of the different single terms. 
		//testingConvectionEqualities();
		
	} else {
		Info << "calculate mean " << endl;
		if (ZoneName != "") { 
			_work = false; 
		} else {
			Info << "Global zone -> calculate " << endl;
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
		} else { 
			tagU_ls = tagU; 
			tagU = U-(*mean_U); 
			tag_p_rgh = p_rgh-(*mean_p_rgh);
			tag_T     = T-(*mean_T);
//			tag_phi  = phi-(*mean_phi);


			update_reynolds();
			update_energy_dUdt();
		}
	}  //.. work
}


void EnergyBalanceTerms::finalize() { 
	if (_work) { 
		if (_calcMean) { 
			finalize_calculate_mean();

			runTime.setTime(mesh.time().startTime(),0); 

			mean_U->write();
			mean_p_rgh->write();
			mean_phi->write();
			mean_T->write(); 

		} else { 
			finalize_calculate_reynolds();

			finalize_energy_balance(); 
		}
	}
}


vector EnergyBalanceTerms::momentum_dUdt() {
	return Integrate(fvc::ddt(U));
}


void EnergyBalanceTerms::finalize_energy_balance() { 
	dimensionedScalar timeSpan = runTime.endTime() - runTime.startTime();
	
	perb_energy_dUdt  /= timeSpan.value(); 
	total_energy_dUdt /= timeSpan.value(); 

	perb_energy_pressure  /= timeSpan.value(); 
	total_energy_pressure /= timeSpan.value(); 

	Info << "Test - Energy - Temporal " << endl;
	Info << "================= " << endl;
	Info << "\t" << total_energy_dUdt << " = " << mean_energy_dUdt << " + " << perb_energy_dUdt << endl;

	Info << "Test - Pressure - Tempora " << endl; 
	Info << "\t" << total_energy_pressure << " = " << mean_energy_pressure << " + " << perb_energy_pressure << endl;



	//Info << "Test - Energy - Mean Balance " << endl;
	//Info << "================= " << endl;

	//volVectorField& meanU 		= *mean_U;
	//surfaceScalarField& meanphi 	= *mean_phi; 
	//volScalarField barEk = 0.5*meanU&meanU;

	// Mean energy convection. 
	//tensor sqrdface_termwise = Integrate( fvc::grad(meanU*barEk) );
	//Info << "\t The mean energy convection:  (total  " << Integrate(meanU&fvc::div(meanphi,meanU))  << ")" << endl; 
	//Info << "\t\t x-flux: meanU*meanEk " << sqrdface_termwise.xx() << endl;
	//Info << "\t\t y-flux: meanU*meanEk " << sqrdface_termwise.yy() << endl; 
	//Info << "\t\t z-flux: meanU*meanEk " << sqrdface_termwise.zz() << endl;  

	// Mean reynolds convection. 
	//Info << "\t The mean reynold flux " << Integrate(fvc::div( reynoldsU& meanU )) << endl; 

}


	// ====================================== dt ====================	
void EnergyBalanceTerms::update_energy_dUdt() { 

	scalar dt = runTime.deltaTValue();
	scalar total_balance = Integrate( U&fvc::ddt(U) )*dt;

	if (runTime.time() > runTime.startTime() ) {

		volVectorField tagU_ddt = (tagU-tagU_ls)/dt;
		perb_energy_dUdt  += Integrate( tagU&tagU_ddt )*dt;
		total_energy_dUdt += total_balance;
	}
}

	// ====================================== pressure ====================	
//- Calculates the <U&grad(p_rgh)> and accumulates it for the average. 
void EnergyBalanceTerms::update_energy_Pressure() {
	scalar dt = mesh.time().deltaTValue();

	total_energy_pressure += Integrate( U&fvc::grad(p_rgh) )*dt;
	perb_energy_pressure  += Integrate( tagU&fvc::grad(tag_p_rgh) )*dt;
} 



	// ====================================== Convection ====================	
//- Calculates the <div(phi,U)>
vector EnergyBalanceTerms::momentum_Convection() {
	return Integrate(fvc::div(phi,U)); 
}

//- Calculates the <U&div(phi,U)> and accumulates it for the average. 
scalar EnergyBalanceTerms::energy_Convection() {
	
	//scalar dt = mesh.time().deltaTValue();

	//scalar local_balance = Integrate( U&fvc::div(phi,U) )*dt;
	//mean_energy_Convection += local_balance;
	return 0;
}

	// ====================================== Diffusion =====================
//- Calculates the <laplacian(AnisitropicDiffusion,U)>
vector EnergyBalanceTerms::momentum_Diffusion() {
	return Integrate(fvc::laplacian(AnisotropicDiffusion,U)); 
}

//- Calculates the <U&laplacian(AnisitropicDiffusion,U)> and accumulates it for the average. 
scalar EnergyBalanceTerms::energy_Diffusion() {

//	scalar dt = mesh.time().deltaTValue();

//	scalar local_balance = Integrate( U&fvc::laplacian(AnisotropicDiffusion,U) )*dt;
//	mean_energy_Diffusion += local_balance;

	return 0;
}

	// ====================================== Pressure  =====================
//- Calculates the <grad(p_rgh)>
vector EnergyBalanceTerms::momentum_Pressure() {
	return Integrate(fvc::grad(p_rgh)); 
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

void  EnergyBalanceTerms::update_reynolds() {
	scalar dt = mesh.time().deltaTValue();
	reynoldsU += tagU*tagU*dt;
}

void EnergyBalanceTerms::finalize_calculate_reynolds() { 
	scalar TotalTime = (runTime.endTime()-runTime.startTime()).value();
	reynoldsU /= TotalTime;
}

void EnergyBalanceTerms::finalize_calculate_mean() {

	scalar TotalTime = (runTime.endTime()-runTime.startTime()).value();

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



// ==================================================================================================================================
// ==================================================================================================================================
// ==================================================================================================================================
// ==================================================================================================================================


//- Checks the equalities in the openFoam energy document (21.1.1)
void EnergyBalanceTerms::testingConvectionEqualities() { 
	// Checks if:  \bar{U}\&fvc::div(\bar{phi},\bar{U})> = <fvc::div(\bar{phi},0.5*\bar{U}\&\bar{U})>
	volVectorField& meanU 		= *mean_U;
	surfaceScalarField& meanphi 	= *mean_phi; 

	volScalarField barEk = 0.5*meanU&meanU;

	tensor sqrdface_termwise = Integrate( fvc::grad(meanU*barEk) );
	Info << "\t"  << sqrdface_termwise.xx() << "," << sqrdface_termwise.yy() << "," << sqrdface_termwise.zz() << " || " << Integrate(meanU&fvc::div(meanphi,meanU)) << endl;

	//volTensorField zeroTensor(diffusionSource,barEk.dimensions());
	tensor diffusionSource; 
    	volTensorField gradU = fvc::grad(meanU); 
    	for(int i=0;i<9;i++) { 
        	diffusionSource.component(i) = Integrate(AnisotropicDiffusion.component(0)*gradU.component(i)*gradU.component(i));
    	}	

	Info << Integrate(meanU&(-fvc::laplacian(AnisotropicDiffusion,meanU))) << " || " << Integrate(AnisotropicDiffusion.component(0)*gradU&&gradU) << endl;
	


	// =========================================================================================================
	//Info << " The <meanU&div(menaphi,meanU)> = <div(meanphi,0.5*meanU&meanU)> equality" << endl;
	//scalar orig = Integrate(meanU&fvc::div(meanphi,meanU));
	//scalar sqredterm = Integrate(fvc::div(meanphi,0.5*meanU&meanU));
	//Info << "\t"  << orig << " - " << sqredterm << " = " << orig-sqredterm << " frac " << (orig-sqredterm)/orig << endl;

	//Info << " The <meanU&div(menaphi,meanU)> = <meanU&grad(meanphi*interpolate(meanU))> equality" << endl;
	//scalar faceterm = Integrate(meanU&fvc::div(meanphi*fvc::interpolate(meanU)));
	//Info << "\t"  << orig << " - " << faceterm << " = " << orig-faceterm << " frac " << (orig-faceterm)/orig << endl;

	//Info << " The <meanU&div(menaphi,meanU)> = <div(meanphi*0.5*interpolate(meanU&meanU))> equality" << endl;
	//scalar sqrdfaceterm = Integrate(fvc::div(0.5*meanphi*fvc::interpolate(meanU&meanU)));
	//Info << "\t"  << orig << " - " << sqrdfaceterm << " = " << orig-sqrdfaceterm << " frac " << (orig-sqrdfaceterm)/orig << endl;
	// =========================================================================================================
}

