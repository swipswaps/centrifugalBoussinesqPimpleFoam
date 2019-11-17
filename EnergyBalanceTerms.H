#ifndef EnergyCalculations_H
#define EnergyCalculations_H 

#include "fvCFD.H"
#include "timeSelector.H"
#include "fvc.H"
#include "fvCFD.H"


#include "meshSearch.H"

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

	volScalarField zoneSelector;
	bool calculateZone; 

	// ========================================================== Averages. 
	volVectorField* mean_U;
	volScalarField* mean_p_rgh;
	volScalarField* mean_T;
	surfaceScalarField* mean_phi;

	scalar	       mean_energy_dUdt; 
	scalar	       mean_energy_Convection;
 	scalar	       mean_energy_Diffusion; 


	// ==========================================================
	// ========================================================== Methods
	// ==========================================================

	void update_mean();
	void finalize_calculate_mean();  


	// ====================================== Temporal ======================
		//- Calculates <ddt(U)>
		inline vector momentum_dUdt(); 

		//- Calculates <U&ddt(U)>
		scalar energy_dUdt(); 

	// ====================================== Convection ====================	
		//- Calculates the <div(phi,U)>
		inline vector momentum_Convection(); 

		//- Calculates the <U&div(phi,U)> and accumulates it for the average. 
		scalar energy_Convection(); 

	// ====================================== Diffusion =====================
		//- Calculates the <laplacian(AnisitropicDiffusion,U)>
		inline vector momentum_Diffusion(); 

		//- Calculates the <U&laplacian(AnisitropicDiffusion,U)> and accumulates it for the average. 
		scalar energy_Diffusion(); 

	// ====================================== Pressure  =====================
		//- Calculates the <grad(p_rgh)>
		inline vector momentum_Pressure(); 

		//- Calculates the <U&grad(p_rgh)> and accumulates it for the average. 
		scalar energy_Pressure(); 

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

};
#endif
