#include "EnergyBalanceTerms.H" 


#define zeroScalar(variable,dims) \
	variable(		\
            IOobject		\
            (			\
                    #variable,		\
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
                    #variable,		\
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
                    #variable,		\
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
                    #variable,		\
                    mesh.time().timeName(),	\
                    mesh,			\
                    IOobject::NO_READ,		\
                    IOobject::AUTO_WRITE	\
            ),					\
            mesh,				\
            dimensionedTensor("zero",dims,tensor::zero)\
	)

#define zeroEnergyScalarTerms(variable,dims,ctype) \
	zeroScalar(mean_##ctype##_##variable,dims),\
	zeroScalar(perb_##ctype##_##variable,dims),\
	zeroScalar(total_##ctype##_##variable,dims)

#define zeroEnergyVectorTerms(variable,dims,ctype) \
	zeroVector(mean_##ctype##_##variable,dims),\
	zeroVector(perb_##ctype##_##variable,dims),\
	zeroVector(total_##ctype##_##variable,dims)


#define zeroEnergyTensorTerms(variable,dims,ctype) \
	zeroTensor(mean_##ctype##_##variable,dims),\
	zeroTensor(perb_##ctype##_##variable,dims),\
	zeroTensor(total_##ctype##_##variable,dims)

#define integrate_variables(term,ctype)\
	scalar	       mean_##ctype##_##term; \
	scalar	       perb_##ctype##_##term;\
	scalar	       total_##ctype##_##term;

#define variables_finalize(term,timespan,ctype)\
	perb_##ctype##_##term  /=timespan;\
	total_##ctype##_##term /=timespan;


/**
	Initializes the internal variables and checks if the mean exists. 
	If the mean exists the operates on a perturbaion calculation mode. 
	Otherwise, calculate the mean. 

*/ 
EnergyBalanceTerms::EnergyBalanceTerms(
				fvMesh& mesh,
				Time&   runTime,
				volVectorField &U,
				surfaceScalarField &phi,
				volTensorField &AnisotropicDiffusion,
				volScalarField &p_rgh,
				volScalarField &T,
				dimensionedScalar beta,
				word ZoneName) :
	_calcMean(false),
	_work(true) ,
	mesh(mesh) ,
	runTime(runTime), 
	U(U) , 
	phi(phi), 
	AnisotropicDiffusion(AnisotropicDiffusion), 
	p_rgh(p_rgh) , 
	T(T)  ,
	U_background(
		IOobject
		(
		    "Umean",
		    runTime.timeName(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::AUTO_WRITE
		),
		mesh
	),
	g(
		IOobject
		(
		    "g",
		    runTime.constant(),
		    mesh,
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	), 
	beta(beta), 
	alpha(dimensionedScalar("alpha",dimless/dimTime,1e-3)),
	zeroScalar(zoneSelector,dimless) ,
	zeroVector(dmeanUdt,U.dimensions()),
	zeroVector(tagU,U.dimensions()),
	zeroVector(tagU_ls,U.dimensions()),
	zeroScalar(tag_p_rgh,p_rgh.dimensions()),
	zeroScalar(tag_T,T.dimensions()),
	zeroScalar(tag_phi,phi.dimensions()) , 
	CenterLookup(mesh.C().internalField().size()) ,
	zeroEnergyScalarTerms(dUdt,dimVelocity*dimVelocity/dimTime,energy),
	zeroEnergyTensorTerms(diffusion,dimVelocity*dimVelocity/dimTime,energy),
	zeroEnergyTensorTerms(gradUsqr,dimVelocity*dimVelocity/dimTime,energy),
	zeroEnergyScalarTerms(pressure,dimVelocity*dimVelocity/dimTime,energy),
	zeroEnergyScalarTerms(nudging,dimVelocity*dimVelocity/dimTime,energy),
	zeroEnergyScalarTerms(potential,dimVelocity*dimVelocity/dimTime,energy),
	zeroTensor(energy_fullFlux, dimVelocity*dimVelocity/dimTime),
	zeroTensor(mean_meanMeanFlux, dimVelocity*dimVelocity/dimTime),
	zeroTensor(mean_perbPerbFlux, dimVelocity*dimVelocity/dimTime),
	zeroTensor(perb_perbMeanFlux, dimVelocity*dimVelocity/dimTime),
	zeroTensor(perb_meanPerbFlux, dimVelocity*dimVelocity/dimTime),
	zeroTensor(perb_perbPerbFlux, dimVelocity*dimVelocity/dimTime),
	zeroTensor(mean_perb_conversion, dimVelocity*dimVelocity/dimTime),
	zeroEnergyTensorTerms(diffusion,dimVelocity/dimTime,momentum),
	zeroEnergyVectorTerms(diffusion,dimVelocity/dimTime,eqn),
	zeroEnergyVectorTerms(diffusion,dimVelocity*dimVelocity/dimTime,eqn_energy),
	zeroTensor(mean_KgradUsqr, dimVelocity*dimVelocity/dimTime),
	zeroTensor(perb_KgradUsqr, dimVelocity*dimVelocity/dimTime)
#ifdef TESTS
       ,zeroScalar(energy_fullFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_meanMeanFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_perbPerbFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(perb_meanPerbFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(perb_perbMeanFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(perb_perbPerbFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(perb_meanMeanFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_meanPerbFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_perbMeanFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_fullFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(perb_fullFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_fullTagFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_fullMeanFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_meanFullFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroScalar(mean_perbFullFlux_eqn,dimVelocity*dimVelocity/dimTime),
	zeroVector(meanPerbFlux_eqn,dimVelocity/dimTime)
#endif 
{

	Info << " =-=-=-=-=- Starting energy terms " << endl;
	// ---------------------------------------------- Setting integration zone 
	globalMeshData globalmeshdata(mesh); 
	const globalIndex& globalindexer = globalmeshdata.globalPointNumbering();

	// Specialized code for the nudging. 
	Info << " \t Nudging  " << endl;
    	meshSearch meshSearch(mesh); 
	scalar CenterOfDomain = (max(mesh.C().component(1))-min(mesh.C().component(1))).value()/2;    

	Foam::meshSearch searchEngine(mesh);
	label centercelli = -1;
	forAll(mesh.C().internalField(), celli) { 
		vector position = mesh.C().internalField()[celli];
		position.component(1) = CenterOfDomain;
		centercelli = searchEngine.findCell(position,centercelli); //findCell(position);
		if (centercelli == -1) { 
			Sout << " Cell not found  in processor "  << Pstream::myProcNo() << ": At cell " 
			     << mesh.C().internalField()[celli] << " label " << celli << " looking for position " << position << " Not Found!" << endl;
		}
		CenterLookup[celli] = centercelli; 
	}

	// --- Specialized code for the nudging : end
	Info << " \t Zone definition " << endl;
	if (ZoneName != "") { 
		label cellzoneID  =  mesh.cellZones().findZoneID(ZoneName);
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
			} else { 
				zoneSelector[currentZoneIndx] = mesh.V()[currentZoneIndx];
			}
		}
	} else { 
		forAll(mesh.V(),cellindx) {
			zoneSelector[cellindx] = mesh.V()[cellindx];
		}
	}

	ZoneVolume = sum(zoneSelector).value(); 
	reduce(ZoneVolume,sumOp<scalar>());

	Info << " The zone volume is " << ZoneVolume <<endl;
	// ---------------------------------------------- 
	IOobject check_mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE);
	if (check_mean_U_Header.headerOk()) { 
		Info << "\t=========================== Mean exists: loading " << endl;
		_calcMean = false;
		Info << "\t\t U" << endl;
		IOobject mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		Info << "\t\t p_rgh" << endl;
		IOobject mean_p_rgh_Header("mean_p_rgh" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		Info << "\t\t T" << endl;
		IOobject mean_T_Header("mean_T" ,mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
		Info << "\t\t phi" << endl;
		IOobject mean_Phi_header("mean_phi",mesh.time().timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);

		mean_U 		= new Foam::volVectorField(mean_U_Header,mesh); 
		mean_p_rgh  = new Foam::volScalarField(mean_p_rgh_Header,mesh);
		mean_T  	= new Foam::volScalarField(mean_T_Header,mesh);
		mean_phi	= new Foam::surfaceScalarField(mean_Phi_header,mesh); 

		Info << "\t\t -- calculating the first u' " << endl;
		tagU = U-(*mean_U); 

		Info << "\t\t Loading the last time step to calculate d[ubar]/dt " << endl;
		word lastTime(name(runTime.endTime().value()));
		dimensionedScalar timeSpan = runTime.endTime() - runTime.startTime();

		Info << "\t\t The timespan is " << timeSpan << endl;

		// get the dubar/dt = (u|1-u|0)/Time
		Info << "\t\t Reading time step " << lastTime;
		volVectorField Ulast(IOobject("U",lastTime,mesh,IOobject::MUST_READ,IOobject::NO_WRITE),
				      mesh);
		Info << " ... Done " << endl; 

		// calculate the mean fields. 
		volVectorField& meanU 		= *mean_U;
		surfaceScalarField& meanphi     = *mean_phi;
		volScalarField& meanT           = *mean_T;

		Info << "\t\t Setting Uc" << endl;
		volVectorField meanUc = meanU;
		setUc(meanUc);

		Info << "\t\t Calculating mean energy terms " << endl;
		Info << "\t\t\t dUdt " << endl;
		mean_energy_dUdt 	= (meanU) & ((Ulast-U)/timeSpan);

		Info << "\t\t\t Pressure " << endl;
		mean_energy_pressure 	= (meanU) & fvc::grad(*mean_p_rgh);

		Info << "\t\t\t nudging " << endl;
		mean_energy_nudging	= alpha*meanU &(U_background-meanUc);

		Info << "\t\t\t potential " << endl;
		mean_energy_potential  = -beta*(g&meanU)()*meanT;  // remember the g is negative!!. 

		Info << "\t\t\t mean mean mean flux " << endl;		
		mean_meanMeanFlux 	  = calculateEnergyFlux(meanphi, meanU, meanU)();
	
		Info << "\t\t\t grad mean U " << endl;		
		volTensorField grad_meanU = fvc::grad(meanU); 

		Info << "\t\t\t (grad mean U)^2 " << endl;		
		forAll(mesh.C(),celli) {
			for (int c=0;c<9;c++) { 
				mean_KgradUsqr[celli].component(c) = AnisotropicDiffusion[celli].component(0)*grad_meanU[celli].component(c)*grad_meanU[celli].component(c);

			}
		}

		Info << "\t\t\t diffusion " << endl;		
		mean_momentum_diffusion 	= calculategradKgrad(meanU)();
		mean_eqn_diffusion 		= fvc::laplacian(AnisotropicDiffusion,meanU);


		label component = 0;
		for(label i=0;i<3;i++) { 
			for (label j=0;j<3;j++,component++) { 
				forAll(mesh.C(),celli) { 
					scalar Ucomponent = meanU[celli].component(component);
					scalar& mean_energy_component = mean_energy_diffusion[celli].component(component); 
					mean_energy_component = Ucomponent*mean_momentum_diffusion[celli].component(component); 
				} //..for cells. 
			}
			forAll(mesh.C(),celli) { 
				scalar& mean_eqn_energy_diffusion_component = mean_eqn_energy_diffusion[celli].component(i);
				mean_eqn_energy_diffusion_component = meanU[celli].component(i)*mean_eqn_diffusion[celli].component(i);
			}
		}

		Info << Integrate(mean_energy_diffusion) << endl;
		Info << "\t\t\t end " << endl;		
		
	} else {
		Info << "calculate mean on global zone" << endl;
		Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;
		_calcMean = true;

		IOobject mean_U_Header("mean_U" ,mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);
		IOobject mean_p_rgh_Header("mean_p_rgh" ,mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);
		IOobject mean_T_Header("mean_T" ,mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);
		IOobject mean_Phi_Header("mean_phi",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE);

		scalar dt = mesh.time().deltaTValue();
		mean_U 		= new Foam::volVectorField(mean_U_Header,U*dt) ; //mesh,dimensionedVector("uu",U.dimensions(),vector::zero)); 
		mean_p_rgh  	= new Foam::volScalarField(mean_p_rgh_Header,mesh,dimensionedScalar("uu",p_rgh.dimensions(),0)); ;
		mean_T  	= new Foam::volScalarField(mean_T_Header,mesh,dimensionedScalar("uu",T.dimensions(),0)); ;
		mean_phi	= new Foam::surfaceScalarField(mean_Phi_Header,mesh,dimensionedScalar("uu",phi.dimensions(),0)); 
	}



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

			tag_phi   = phi-*mean_phi;

			Info << " Updating dudt " << endl;			
			update_energy_dUdt();

			Info << " Updating pressure " << endl;
			update_energy_Pressure();

			Info << " Updating convection " << endl;
			update_energy_Convection();

			Info << " Updating nudging " << endl;
			update_energy_Nudging();  

			Info << " Updating potential " << endl;
			update_energy_Potential();

			Info << " Updating diffusion " << endl;
			update_energy_Diffusion();

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
			dimensionedScalar timeSpan = runTime.endTime() - runTime.startTime();
			finalize_calculate_perb(); 

			variables_finalize(diffusion,timeSpan.value(),momentum);
			variables_finalize(diffusion,timeSpan.value(),eqn);
			variables_finalize(diffusion,timeSpan.value(),eqn_energy);

			// Tests
				Info << " Tests " << endl;
				Info << "=================" << endl;
				test_energy_dUdt();
				test_energy_Pressure();
				test_energy_Diffusion();
				test_energy_Convection();
				test_energy_Nudging();
				test_enetgy_Potential();

		}
	}


}


void EnergyBalanceTerms::finalize_calculate_perb() { 

	dimensionedScalar timeSpan = runTime.endTime() - runTime.startTime();
	
	variables_finalize(dUdt,timeSpan.value(),energy); 
	variables_finalize(pressure,timeSpan.value(),energy); 
	variables_finalize(nudging,timeSpan.value(),energy);
	variables_finalize(diffusion,timeSpan.value(),energy);
	variables_finalize(potential,timeSpan.value(),energy);        

	energy_fullFlux		/= timeSpan;
	mean_meanMeanFlux	/= timeSpan;
	mean_perbPerbFlux	/= timeSpan;
	perb_perbMeanFlux	/= timeSpan;
	perb_meanPerbFlux	/= timeSpan;
	perb_perbPerbFlux	/= timeSpan;

	mean_perb_conversion    /= timeSpan; 
	perb_KgradUsqr		/= timeSpan; 

#ifdef TESTS
	energy_fullFlux_eqn    /= timeSpan;   // U&fvc::div(phi,U);
	mean_meanMeanFlux_eqn    /= timeSpan; // meanU&fvc::div(meanphi,meanU);
	mean_perbPerbFlux_eqn    /= timeSpan; // tagU&fvc::div(tag_phi,meanU);
	perb_meanPerbFlux_eqn    /= timeSpan; // tagU&fvc::div(meanphi,tagU);   
	perb_perbMeanFlux_eqn    /= timeSpan; // tagU&fvc::div(tag_phi,meanU);
	perb_perbPerbFlux_eqn    /= timeSpan; // tagU&fvc::div(tag_phi,tagU);    

	mean_fullFlux_eqn    /= timeSpan;
	perb_fullFlux_eqn    /= timeSpan;
	
	mean_fullTagFlux_eqn /= timeSpan;
	mean_fullMeanFlux_eqn/= timeSpan;
	mean_meanFullFlux_eqn/= timeSpan;
	mean_perbFullFlux_eqn/= timeSpan;

	// close to zero.
	perb_meanMeanFlux_eqn    /= timeSpan;
	mean_meanPerbFlux_eqn    /= timeSpan;
	mean_perbMeanFlux_eqn    /= timeSpan;
#endif
}


	// ====================================== dt ====================	
void EnergyBalanceTerms::update_energy_dUdt() { 

	dimensionedScalar dim_dt = runTime.deltaT();
	scalar dt = runTime.deltaTValue();
	if (runTime.time() > runTime.startTime() ) {

		volVectorField tagU_ddt = (tagU-tagU_ls)/dim_dt;
		perb_energy_dUdt  += tagU&tagU_ddt*dt;
		total_energy_dUdt += ( U&fvc::ddt(U) )*dt;
	}
}

void  EnergyBalanceTerms::test_energy_dUdt() { 

	scalar mean  = Integrate(mean_energy_dUdt); 
	scalar perb  = Integrate(perb_energy_dUdt); 
	scalar total = Integrate(total_energy_dUdt); 

 	Info << " u*du/dt " << total << " = " << mean << " + " << perb << " = " << perb+mean << " Relative " << (perb+mean)/total << endl;
 	
 	Info << " =========perb========== " << endl; 
	forAll(perb_energy_dUdt,celli) { 
			Info << celli << " "  << perb_energy_dUdt[celli] << endl; 
			
			if (celli > 10) break;
	}
 	Info << " =================== " << endl; 
 	Info << " =========mean========== " << endl; 
	forAll(mean_energy_dUdt,celli) { 
			Info << celli << " "  << mean_energy_dUdt[celli] << endl; 
			
			if (celli > 10) break;
	}
 	Info << " =================== " << endl; 
 	Info << " =========total========== " << endl; 
	forAll(total_energy_dUdt,celli) { 
			Info << celli << " "  << total_energy_dUdt[celli] << endl; 
			
			if (celli > 10) break;
	}
 	Info << " =================== " << endl; 
 	
}

	// ====================================== pressure ====================	

//- Calculates the <U&grad(p_rgh)> and accumulates it for the average. 
void EnergyBalanceTerms::update_energy_Pressure() {
	scalar dt = mesh.time().deltaTValue();

	total_energy_pressure +=  U&fvc::grad(p_rgh)*dt;
	volVectorField grad_tag_p_rgh = fvc::grad(p_rgh)-fvc::grad(*mean_p_rgh);
	perb_energy_pressure  +=  tagU&grad_tag_p_rgh*dt;
} 

void  EnergyBalanceTerms::test_energy_Pressure() { 

	scalar mean  = Integrate(mean_energy_pressure); 
	scalar perb  = Integrate(perb_energy_pressure); 
	scalar total = Integrate(total_energy_pressure); 

 	Info << " u*pressure " << total << " = " << mean << " + " << perb << " = " << perb+mean << " Relative " << (perb+mean)/total << endl;
}



	// ====================================== diffusion ====================	
void  EnergyBalanceTerms::update_energy_Diffusion() {


	scalar dt = mesh.time().deltaTValue();

	volTensorField total_momentum_current = calculategradKgrad(U)();

	volVectorField total_eqn_diffusion_current = fvc::laplacian(AnisotropicDiffusion,U);
	volVectorField perb_eqn_diffusion_current = fvc::laplacian(AnisotropicDiffusion,tagU);

	total_momentum_diffusion       += total_momentum_current*dt; 
	total_eqn_diffusion 	       += total_eqn_diffusion_current*dt;

	volTensorField perb_momentum_current = calculategradKgrad(tagU)(); 
	perb_momentum_diffusion        += perb_momentum_current*dt; 
	perb_eqn_diffusion 	       += perb_eqn_diffusion_current*dt;


	label component = 0;
	for(label i=0;i<3;i++) { 
		for (label j=0;j<3;j++,component++) { 
			forAll(mesh.C(),celli) { 
				scalar Ucomponent = U[celli].component(j);
				scalar tagUcomponent = tagU[celli].component(j);

				scalar& total_energy_component 	 = total_energy_diffusion[celli].component(component); 
				total_energy_component 		+= Ucomponent*total_momentum_current[celli].component(component)*dt; 

				scalar& perb_energy_component    = perb_energy_diffusion[celli].component(component); 
				perb_energy_component 		+= tagUcomponent*perb_momentum_current[celli].component(component)*dt; 
			} //..for cells. 

		}


		forAll(mesh.C(),celli) { 
			scalar& total_eqn_energy_diffusion_component = total_eqn_energy_diffusion[celli].component(i);
			total_eqn_energy_diffusion_component += U[celli].component(i)*total_eqn_diffusion_current[celli].component(i)*dt;

			scalar& perb_eqn_energy_diffusion_component = perb_eqn_energy_diffusion[celli].component(i);
			perb_eqn_energy_diffusion_component += tagU[celli].component(i)*perb_eqn_diffusion_current[celli].component(i)*dt;

		}

	}

	Info << "\t\t\t (grad tagU)^2 " << endl;		
	volTensorField grad_perbU = fvc::grad(tagU); 

	forAll(mesh.C(),celli) {
		for (int c=0;c<9;c++) { 
			perb_KgradUsqr[celli].component(c) += AnisotropicDiffusion[celli].component(0)*grad_perbU[celli].component(c)*grad_perbU[celli].component(c)*dt;

		}
	}


	// calculating the (grad U)^2. 
//	perb_energy_gradUsqr  += KgradtagUsqr*dt;
//	total_energy_gradUsqr += KgradUsqr*dt;  
}

void  EnergyBalanceTerms::test_energy_Diffusion() { 

	tensor mean  = Integrate(mean_energy_diffusion); 
	tensor perb  = Integrate(perb_energy_diffusion); 
	tensor total = Integrate(total_energy_diffusion); 

	tensor momentum_mean  = Integrate(mean_momentum_diffusion); 
	tensor momentum_perb  = Integrate(perb_momentum_diffusion); 
	tensor momentum_total = Integrate(total_momentum_diffusion); 

	vector eqn_total      = Integrate(total_eqn_diffusion); 
	vector eqn_mean       = Integrate(mean_eqn_diffusion); 
	vector eqn_perb       = Integrate(perb_eqn_diffusion); 

	vector eqn_energy_total      = Integrate(total_eqn_energy_diffusion); 
	vector eqn_energy_mean       = Integrate(mean_eqn_energy_diffusion); 
	vector eqn_energy_perb       = Integrate(perb_eqn_energy_diffusion); 

	word component;
	word dir1;
	word dir2;

 	Info << " u*diffusion " << endl;
	Info << " \t---- momentum Eqn. terms " << endl;
	Info << "\ttotal laplacian(K,u): " << eqn_total.component(0) << " = " << momentum_total.component(0) + momentum_total.component(3) + momentum_total.component(6) << endl;
	Info << "\ttotal laplacian(K,v): " << eqn_total.component(1) << " = " << momentum_total.component(1) + momentum_total.component(4) + momentum_total.component(7) << endl;
	Info << "\ttotal laplacian(K,w): " << eqn_total.component(2) << " = " << momentum_total.component(2) + momentum_total.component(5) + momentum_total.component(8) << endl;
	Info << endl;
	Info << "\tmean laplacian(K,u): " << eqn_mean.component(0) << " = " << momentum_mean.component(0) + momentum_mean.component(3) + momentum_mean.component(6) << endl;
	Info << "\tmean laplacian(K,v): " << eqn_mean.component(1) << " = " << momentum_mean.component(1) + momentum_mean.component(4) + momentum_mean.component(7) << endl;
	Info << "\tmean laplacian(K,w): " << eqn_mean.component(2) << " = " << momentum_mean.component(2) + momentum_mean.component(5) + momentum_mean.component(8) << endl;
	Info << endl;
	Info << "\tperb laplacian(K,u): " << eqn_perb.component(0) << " = " << momentum_perb.component(0) + momentum_perb.component(3) + momentum_perb.component(6) << endl;
	Info << "\tperb laplacian(K,v): " << eqn_perb.component(1) << " = " << momentum_perb.component(1) + momentum_perb.component(4) + momentum_perb.component(7) << endl;
	Info << "\tperb laplacian(K,w): " << eqn_perb.component(2) << " = " << momentum_perb.component(2) + momentum_perb.component(5) + momentum_perb.component(8) << endl;

	Info << " \t---- momentum " << endl;
	label c=0;
	for (int i=0 ; i < 3 ; i++ ) {
		switch (i) {
			case 0: 
				dir1 = "x";
				break;
			case 1: 
				dir1 = "y";
				break;
			case 2: 
				dir1 = "z";
				break;
		};
		for (int j=0 ; j < 3 ; j++,c++ ) {
			switch (j) {
				case 0: 
					component="u";
					dir2 = "x";
					break;
				case 1: 
					component="v";
					dir2 = "y";
					break;
				case 2: 
					component="w";
					dir2 = "z";
					break;
			};



			Info << "\t\t d" << component << "/d"<< dir1 << dir1 << 
					" total " << momentum_total.component(c) << " = " 
						  << momentum_mean.component(c)  << " + " << momentum_perb.component(c) << " = " << momentum_perb.component(c)+momentum_mean.component(c);

					if (fabs(momentum_total.component(c)) >  1e-5) { 
						  Info << " Relative " << (momentum_perb.component(c)+momentum_mean.component(c))/momentum_total.component(c) << endl;
					} else { 
						Info << endl;
					}
		}
	}
	
	Info << " \t---- energy Eqn. terms " << endl;
	Info << "\ttotal laplacian(K,u): " << eqn_energy_total.component(0) << " = " << total.component(0) + total.component(3) + total.component(6) << endl;
	Info << "\ttotal laplacian(K,v): " << eqn_energy_total.component(1) << " = " << total.component(1) + total.component(4) + total.component(7) << endl;
	Info << "\ttotal laplacian(K,w): " << eqn_energy_total.component(2) << " = " << total.component(2) + total.component(5) + total.component(8) << endl;
	Info << endl;
	Info << "\tmean laplacian(K,u): " << eqn_energy_mean.component(0) << " = " << mean.component(0) + mean.component(3) + mean.component(6) << endl;
	Info << "\tmean laplacian(K,v): " << eqn_energy_mean.component(1) << " = " << mean.component(1) + mean.component(4) + mean.component(7) << endl;
	Info << "\tmean laplacian(K,w): " << eqn_energy_mean.component(2) << " = " << mean.component(2) + mean.component(5) + mean.component(8) << endl;
	Info << endl;
	Info << "\tperb laplacian(K,u): " << eqn_energy_perb.component(0) << " = " << perb.component(0) + perb.component(3) + perb.component(6) << endl;
	Info << "\tperb laplacian(K,v): " << eqn_energy_perb.component(1) << " = " << perb.component(1) + perb.component(4) + perb.component(7) << endl;
	Info << "\tperb laplacian(K,w): " << eqn_energy_perb.component(2) << " = " << perb.component(2) + perb.component(5) + perb.component(8) << endl;

	Info << " \t---- energy " << endl;
	c=0;
	for (int i=0 ; i < 3 ; i++ ) {
		switch (i) {
			case 0: 
				dir1 = "x";
				break;
			case 1: 
				dir1 = "y";
				break;
			case 2: 
				dir1 = "z";
				break;
		};
		for (int j=0 ; j < 3 ; j++,c++ ) {
			switch (j) {
				case 0: 
					component="u";
					dir2 = "x";
					break;
				case 1: 
					component="v";
					dir2 = "y";
					break;
				case 2: 
					component="w";
					dir2 = "z";
					break;
			};
			Info << "\t\t d" << component << "/d"<< dir1 << dir1 << 
					" total " << total.component(c) << " = " << mean.component(c) << " + " << perb.component(c) << " = " << perb.component(c)+mean.component(c);
					if (fabs(total.component(c)) >  1e-5) { 
						  Info << " Relative " << (perb.component(c)+mean.component(c))/total.component(c) << endl;
					} else { 
						Info << endl;
					}
		}
	}
}


	// ====================================== Convection ====================	
void EnergyBalanceTerms::update_energy_Convection() {
	
	scalar dt = mesh.time().deltaTValue();
	
	volVectorField& meanU 		= *mean_U;
	surfaceScalarField& meanphi     = *mean_phi;

	// ----------
	// u div(u*u) 
	tmp<volTensorField> energy_fullFlux_current = calculateEnergyFlux(phi, U, U);
	energy_fullFlux   += energy_fullFlux_current()*dt;

	// ------------- 			-----
	// ubar div(u'*u') ==  div[ubar *Etag_k] + u'u'div(ubar) == convection of Reynold by mean + conversion terms. 
	tmp<volTensorField> mean_perbPerbFlux_current = calculateEnergyFlux(tag_phi, tagU, meanU);
	mean_perbPerbFlux += mean_perbPerbFlux_current()*dt;

	//		       --------------   -------------
	// u' div(u'*ubar) ==  u'ubar div[u'] + u'u'div(ubar) == Conversion terms. 
	// !!! ----- Note: Calculates the sum of the conversion terms per wind component ---!!!
	tmp<volTensorField> conversionMethodII = calculateEnergyFlux(tag_phi, meanU, tagU);
	perb_perbMeanFlux += conversionMethodII()*dt;

	// -----------------------
	// u'_j d(ubar'_i u'_j)/dxi == mean convection of perb Kinetic energy. 
	tmp<volTensorField> perb_meanPerbFlux_current = calculateEnergyFlux(meanphi, tagU, tagU);
	perb_meanPerbFlux += perb_meanPerbFlux_current()*dt;

	// ----------
	// u' div(u'*u') == perturbation convection of perb kinetic energy 
	tmp<volTensorField> perb_perbPerbFlux_current = calculateEnergyFlux(tag_phi, tagU, tagU);
	perb_perbPerbFlux += perb_perbPerbFlux_current()*dt;

	/// Calculating conversion using method I, calculates all the conversion terms. 
	volTensorField reynolds  = tagU*tagU; 
	volTensorField grad_MeanU = fvc::grad(meanU); 
	
	for(int i=0; i<9; i++) { 
		forAll(mesh.C(),celli) { 
			scalar& mean_perb_conversion_component = mean_perb_conversion[celli].component(i);
			mean_perb_conversion_component        += reynolds[celli].component(i)*grad_MeanU[celli].component(i)*dt;
		}

	}



#ifdef TESTS
//	test_energy_Convection_Mean_perbperb(mean_perbPerbFlux_current); 
//	test_energy_Convection_perb_perbMeanFlux(conversionMethodII);
//	test_energy_Convection_perb_meanPerbFlux(perb_meanPerbFlux_current);
//	test_energy_Convection_perb_perbPerbFlux(perb_perbPerbFlux_current);
	update_energy_Convection_eqn_sum();
#endif

}

void   EnergyBalanceTerms::test_energy_Convection() { 

	word component;
	word dir1;
	word dir2;
	
	Info << "Conversion terms " << endl;

	tensor ConversionTerms_methodI  = Integrate(perb_perbMeanFlux); /// method I only compares sum!. 
	tensor ConversionTerms_methodII = Integrate(mean_perb_conversion);  

	scalar tagUdiv                  = Integrate(fvc::div(tag_phi)); 
	Info << "\t note that the error in the following terms can be of order tag U " << tagUdiv << endl;
	Info << "\tubar " << ConversionTerms_methodI.component(0)  + ConversionTerms_methodI.component(3)  + ConversionTerms_methodI.component(6) << "=" 
			<< ConversionTerms_methodII.component(0) + ConversionTerms_methodII.component(3) + ConversionTerms_methodII.component(6)  << endl;
	Info << "\tvbar " << ConversionTerms_methodI.component(1)  + ConversionTerms_methodI.component(4)  + ConversionTerms_methodI.component(7) << "=" 
			<< ConversionTerms_methodII.component(1) + ConversionTerms_methodII.component(4) + ConversionTerms_methodII.component(7)  << endl;
	Info << "\twbar " << ConversionTerms_methodI.component(2)  + ConversionTerms_methodI.component(5)  + ConversionTerms_methodI.component(8) << "=" 
			<< ConversionTerms_methodII.component(2) + ConversionTerms_methodII.component(5) + ConversionTerms_methodII.component(8)  << endl;

	Info << "Direct calculation of each term " << endl; 
	Info <<" u'u'dU/dx " << ConversionTerms_methodII.component(0) << endl;
	Info <<" u'v'dU/dy " << ConversionTerms_methodII.component(3) << endl;
	Info <<" u'w'dU/dz " << ConversionTerms_methodII.component(6) << endl;

	Info <<" v'u'dV/dx " << ConversionTerms_methodII.component(1) << endl;
	Info <<" v'v'dV/dy " << ConversionTerms_methodII.component(4) << endl;
	Info <<" v'w'dV/dz " << ConversionTerms_methodII.component(7) << endl;

	Info <<" w'u'dV/dx " << ConversionTerms_methodII.component(2) << endl;
	Info <<" w'v'dV/dy " << ConversionTerms_methodII.component(5) << endl;
	Info <<" w'w'dV/dz " << ConversionTerms_methodII.component(8) << endl;
	

}


	// ====================================== Nudging ====================	

void EnergyBalanceTerms::update_energy_Nudging() { 
	scalar dt = mesh.time().deltaTValue();
	volVectorField tagUc = tagU;
	setUc(tagUc); 
	volVectorField Uc = U;
	setUc(Uc); 

	perb_energy_nudging  += -(alpha*tagU&tagUc)*dt;
	total_energy_nudging +=  alpha*U&(U_background-Uc)*dt;

}


void EnergyBalanceTerms::test_energy_Nudging() { 

	scalar total_nudging = Integrate(total_energy_nudging); 
	scalar perb_nudging  = Integrate(perb_energy_nudging); 
	scalar mean_nudging  = Integrate(mean_energy_nudging); 

	Info << "Nudging " << total_nudging  << " = " << mean_nudging+perb_nudging << " = " << mean_nudging << "+ " <<  perb_nudging << endl;
 
}

	// ====================================== Potential ====================	

void EnergyBalanceTerms::update_energy_Potential() {   //- Calculates the potential energy 
	scalar dt = mesh.time().deltaTValue();
	// remember the g is negative!!. 
	perb_energy_potential 	-= beta*tag_T*(g&tagU)()*dt; 
	total_energy_potential  -= beta*T*(g&U)()*dt;

}

void EnergyBalanceTerms::test_enetgy_Potential() { 

	scalar total_potential = Integrate(total_energy_potential); 
	scalar perb_potential  = Integrate(perb_energy_potential); 
	scalar mean_potential  = Integrate(mean_energy_potential); 

	Info << "Potential " << total_potential  << " = " << mean_potential+perb_potential << " = " << mean_potential << "+ " <<  perb_potential << endl;
}

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

	scalar dt = mesh.time().deltaTValue();
	scalar TotalTime = (runTime.endTime()-runTime.startTime()).value()+dt;

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
tmp<volTensorField> EnergyBalanceTerms::calculateFlux(surfaceScalarField& iphi, volVectorField& iU) { 
	
	surfaceVectorField interpolateU		=	fvc::interpolate(iU);
	surfaceTensorField interpolateU2	= 	(iphi*(mesh.Sf()/mesh.magSf()))*interpolateU;

	// correction of the boundary 
	forAll(mesh.boundary(), patchi)
	{	
		Foam::fvsPatchField<Foam::tensor >& boundary_interpolateU2 	= interpolateU2.boundaryField()[patchi];

		const Foam::fvsPatchField<scalar>& boundary_phi	     = iphi.boundaryField()[patchi];
		const Foam::fvsPatchField<Foam::vector>& boundary_interpolateU    = interpolateU.boundaryField()[patchi];
		const Foam::fvPatch& boundary_mesh 	     = mesh.boundary()[patchi];

		forAll(boundary_mesh, facei)
        	{
		    vector n_inside = (boundary_mesh.Sf()[facei]/boundary_mesh.magSf()[facei]);
		    for (int i=0;i<3;i++) { 
			scalar val = n_inside.component(i);
			n_inside.component(i) = fabs(val); 
		    }
	
		    tensor tt = (boundary_phi[facei]*n_inside )*boundary_interpolateU[facei]; 
		    for (int i=0; i<9;i++) {  
				boundary_interpolateU2[facei].component(i) = tt.component(i);
		    }
	        }

 	 }
	return fvc::surfaceIntegrate(interpolateU2);
}

/*	The flux are:  [_{x,y,z} is the derivative in that direction].
//
//			/ (uu)_x  (uv)_x (uw)_x \
//			| (uv)_y  (vv)_y (vw)_y |
//			\ (uw)_z  (vw)_z (ww)_z /	 
//
//	The velocity components are IUenergy = [a,b,c]
//	The output is : 
//			/ a*(uu)_x  b*(uv)_x c*(uw)_x \
//			| a*(uv)_y  b*(vv)_y c*(vw)_y |
//			\ a*(uw)_z  b*(vw)_z c*(ww)_z /	 
*/			
tmp<volTensorField> EnergyBalanceTerms::calculateEnergyFlux(tmp<volTensorField> iFluxTensorPtr, volVectorField& iUEnergy) {

	tmp<volTensorField> retPtr(new volTensorField(
					    IOobject		
					    (			
						    "tmpflux",		
						    mesh.time().timeName(),	
						    mesh,			
						    IOobject::NO_READ,		
						    IOobject::AUTO_WRITE	
					    ),					
					    mesh,				
					    dimensionedTensor("zero",dimVelocity*dimVelocity/dimTime,tensor::zero)
					)
				); 

	volTensorField& iFluxTensor = iFluxTensorPtr();
	volTensorField& ret = retPtr();
	int component = 0; 
	for (int i=0; i<3;i++) { 
		for (int j=0; j<3;j++,component++) {  
			forAll(mesh.C(),celli) { 
				scalar& ret_component = ret[celli].component(component);
			 	ret_component = iFluxTensor[celli].component(component) * iUEnergy[celli].component(j);
			}
		}
	}

	return retPtr;
}

//   calculates iUEnergy*div(phi,u)  
//	 so the energy terms is 
//      <iUEnergy type>_<phi type><u type> 
//    for example: 
// 			calculateEnergyFlux(tag_phi, mean_U, tag_u) is actually perb_perbmean 
//
tmp<volTensorField> EnergyBalanceTerms::calculateEnergyFlux(surfaceScalarField& iphi, volVectorField& iU, volVectorField& iUEnergy) { 
	tmp<volTensorField> ret =  calculateEnergyFlux(  
					calculateFlux(iphi,iU),
					iUEnergy
				   );
	return ret; 
}

void EnergyBalanceTerms::setUc(volVectorField& iU) { 
	forAll(mesh.C(), celli) { 
		iU[celli] = iU[CenterLookup[celli]]; 
	}
}

tmp<volTensorField>  EnergyBalanceTerms::calculategradKgrad(volVectorField& iU) { 
	tmp<volTensorField> ret(new volTensorField(
					    IOobject		
					    (			
						    "tmplaplacian",		
						    mesh.time().timeName(),	
						    mesh,			
						    IOobject::NO_READ,		
						    IOobject::AUTO_WRITE	
					    ),					
					    mesh,				
					    dimensionedTensor("zero",dimVelocity/dimTime,tensor::zero)
					)
				); 

	const surfaceTensorField gamma = fvc::interpolate(AnisotropicDiffusion); 
	const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
	surfaceTensorField KsnGradU = (mesh.Sf() & gamma)*fvc::snGrad(iU);
	forAll(mesh.boundary(), patchi)
    	{
		forAll(mesh.boundary()[patchi], facei)
		{
			label component=0;
			for(int i=0; i<3 ; i++) { 
				scalar dir = Sn.boundaryField()[patchi][facei].component(i); 
				for(int j=0; j<3 ; j++,component++) { 
					KsnGradU.boundaryField()[patchi][facei].component(component) *= dir;
				}
			}
		}
	}
	ret() = fvc::div(KsnGradU);	
	return ret;
}

// ==================================================================================================================================
// ==================================================================================================================================



void EnergyBalanceTerms::checkMomentumBalance_Timestep() { 
	
	Info << "\t\tMomentum balance for time step " << runTime.timeName() <<  endl;
	Info << "---------------------------------------------------------------------" << endl; 

	volVectorField Uc = U;
	setUc(Uc);
	volVectorField NudgingTerm = alpha*(U_background-Uc);
	
	vector TotalMomentumIntegration = Integrate(fvc::ddt(U) + fvc::div(phi, U) + NudgingTerm - fvc::laplacian(AnisotropicDiffusion,U) + fvc::grad(p_rgh) + g*beta*T);
	Info << "\t Momentum conservation " << TotalMomentumIntegration<< endl;
}

#ifdef TESTS
#include "EnergyBalanceTests.C"
#endif

