#ifndef EnergyCalculations_H
#define EnergyCalculations_H 

#include "fvCFD.H"
#include "timeSelector.H"
#include "fvc.H"
#include "fvCFD.H"


#include "meshSearch.H"


#define energy_variables(term)\
	scalar	       mean_energy_##term; \
	scalar	       perb_energy_##term;\
	scalar	       total_energy_##term;
	

#define energy_variables_set_zero(term)\
	mean_energy_##term =0; \
	perb_energy_##term =0;\
	total_energy_##term=0;


#define energy_variables_finalize(term,timespan)\
	perb_energy_##term  /=timespan;\
	total_energy_##term /=timespan;


using namespace Foam;

/**
	< > - spatial averaging. 

**/

class EnergyBalanceTerms { 

	public:


	EnergyBalanceTerms(	
				fvMesh& mesh,
				Time&   runTime,
				volVectorField &U,
				surfaceScalarField &phi,
				volTensorField &AnisotropicDiffusion,
				volScalarField &p_rgh,
				volScalarField &T,
				word ZoneName=""
			   );



	void update();
	void finalize();

	protected:

	bool _calcMean;
	bool _work;    // if mean is not present and ZoneName is not "", don't do anything. 


	// ==========================================================
	// ========================================================== Fields.
	fvMesh&	  mesh; 	
	Time  &   runTime; 

	// ==========================================================
	
	volVectorField & U;
	surfaceScalarField &phi;
	volTensorField &AnisotropicDiffusion;
	volScalarField &p_rgh;
	volScalarField &T;

	bool calculateZone; 
	volScalarField zoneSelector;


	// ========================================================== Averages. 
	volVectorField* mean_U;
	volScalarField* mean_p_rgh;
	volScalarField* mean_T;
	surfaceScalarField* mean_phi;
	
	// ------------ mean, perb, total 	
	energy_variables(dUdt); 
	energy_variables(pressure); 
	energy_variables(convection); 


	// ----------- terms termwise - mean 
	tensor mean_convection_termwise;

	volVectorField dmeanUdt;
	volVectorField tagU;
	volVectorField tagU_ls;

	volScalarField tag_p_rgh;
	volScalarField tag_T;
//	surfaceScalarField tag_phi;

	volTensorField reynoldsU; 

	// ==========================================================
	// ========================================================== Methods
	// ==========================================================

	// ------------ calculating the mean stage. 
	void update_mean();
	void finalize_calculate_mean();  	


	// ------------ calculating the perturbations. 
	void update_reynolds();
	void update_energy_dUdt(); //- Calculates <U&ddt(U)> and the perturbation. 
	void update_energy_Pressure(); //- Calculates the <U&grad(p_rgh)> and the perturbation part. 

	void finalize_calculate_reynolds();   // finish calculation of the reynolds numbers. 
	void finalize_energy_balance();


	// ====================================== Temporal ======================

	// ====================================== Convection ====================	

		//- Calculates the <U&div(phi,U)> and accumulates it for the average. 
		scalar energy_Convection(); 

	// ====================================== Diffusion =====================

		//- Calculates the <U&laplacian(AnisitropicDiffusion,U)> and accumulates it for the average. 
		scalar energy_Diffusion(); 

	// ====================================== Pressure  =====================


	// ====================================== Nudging =======================




	// ==========================================================
	// ========================================================== Variables. 
	// ==========================================================


	// ==========================================================
	// ========================================================== util. 
	// ==========================================================

	scalar Integrate(const volScalarField & field) { 
			scalar sfield;
			if(calculateZone){
				sfield = sum(field*zoneSelector).value(); 
			} else { 
				sfield = sum(field).value(); 
			}
			reduce(sfield,sumOp<scalar>());
			return sfield;
	}

	scalar Integrate(const tmp<volScalarField> & field) { 
			scalar sfield;
			if(calculateZone){
				sfield = sum(field*zoneSelector).value(); 
			} else { 
				sfield = sum(field).value(); 
			}
			reduce(sfield,sumOp<scalar>());
			return sfield;
	}

	vector Integrate(const volVectorField & field) { 
			vector sfield;
			if(calculateZone){
				sfield = sum(field*zoneSelector).value(); 
			} else { 
				sfield = sum(field).value(); 
			}
			reduce(sfield,sumOp<vector>());
			return sfield;
	}


	vector Integrate(const tmp<volVectorField> & field) { 
			vector sfield;
			if(calculateZone){
				sfield = sum(field*zoneSelector).value(); 
			} else { 
				sfield = sum(field).value(); 
			}
			reduce(sfield,sumOp<vector>());
			return sfield;
	}

	tensor Integrate(const tmp<volTensorField> & field) { 
			tensor sfield;
			if(calculateZone){
				sfield = sum(field*zoneSelector).value(); 
			} else { 
				sfield = sum(field).value(); 
			}
			reduce(sfield,sumOp<tensor>());
			return sfield;
	}

	tensor Integrate(const volTensorField& field) { 
			tensor sfield;
			if(calculateZone){
				sfield = sum(field*zoneSelector).value(); 
			} else { 
				sfield = sum(field).value(); 
			}
			reduce(sfield,sumOp<tensor>());
			return sfield;
	}

	// =======================================================================

	//- Checks the equalities in the openFoam energy document (21.1.1)
	void testingConvectionEqualities(); 

};
#endif
