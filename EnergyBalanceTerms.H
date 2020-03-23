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


#define energy_fields(type,term)\
	vol##type##Field mean_energy_##term; \
	vol##type##Field perb_energy_##term;\
	vol##type##Field total_energy_##term;
	
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
				dimensionedScalar beta,
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

	volVectorField U_background; // The background wind. In the centerfuginal code it is referred to as Umean. 
	uniformDimensionedVectorField g; 
	dimensionedScalar beta;
	const dimensionedScalar alpha; // the nudging coefficient. 1e-3/s
	
	bool calculateZone; 
	volScalarField zoneSelector;
	double ZoneVolume;


	// ========================================================== Averages. 
	volVectorField* mean_U;
	volScalarField* mean_p_rgh;
	volScalarField* mean_T;
	surfaceScalarField* mean_phi;


	// --- tag Fields 
	volVectorField dmeanUdt;
	volVectorField tagU;
	volVectorField tagU_ls;

	volScalarField tag_p_rgh;
	volScalarField tag_T;
	surfaceScalarField tag_phi;

	List<int> CenterLookup; 

	// - Calculated fields - mean 
	energy_fields(Scalar,dUdt);    		// dt terms. 
	energy_fields(Tensor,diffusion);    	// diffusion terms.
	energy_fields(Tensor,gradUsqr);    	// The diffusion terms after the integration. 
	energy_fields(Scalar,pressure);    	// pressure terms.
	energy_fields(Scalar,nudging);    	// nudging terms.
	energy_fields(Scalar,potential);    	// potential energy terms.        

	// - Convection terms. 
	volTensorField energy_fullFlux; 	// total flux of energy. 
	volTensorField mean_meanMeanFlux; 	// mean flux of mean kinetic energy. 
	volTensorField mean_perbPerbFlux;       // = Conversion term + reynolds flux. 
	volTensorField perb_perbMeanFlux;       // = Conversion term (Method II) . 
	volTensorField perb_meanPerbFlux; 	// = Mean Flux of perturbation kinetic energy. 
	volTensorField perb_perbPerbFlux; 	// = Perturbation flux of perturbation kinetic energy. 

	volTensorField mean_perb_conversion; 	// = Direct calculation of the conversion term. (Method I).


	// ==========================================================
	// ========================================================== Methods
	// ==========================================================

	// ------------ calculating the mean stage. 
	void update_mean();


	// ------------ calculating the perturbations. 
	void update_energy_dUdt(); //- Calculates <U&ddt(U)> and the perturbation. 
	void update_energy_Pressure(); //- Calculates the <U&grad(p_rgh)> and the perturbation part. 
	void update_energy_Nudging();  //- Calculates the nudging force. 
	void update_energy_Potential();  //- Calculates the potential energy 
	void update_energy_Convection();  //- Calculates the total convection
	void update_energy_Diffusion();  //- Calculates the diffusion

	void finalize_calculate_perb();
	void finalize_calculate_mean();  	

	// ==========================================================
	// ========================================================== Variables. 
	// ==========================================================


	// - Calculates the termwise flux with the surface flux iphi and the cell velocity iU. 
	tmp<volTensorField> calculateFlux(surfaceScalarField& iphi, volVectorField& iU) ;

	// - Multiplies the tensor flux by iUEnergy. 
	// - iUEnergy.coponent(0) --> line 0 in the tensor. 
	// - iUEnergy.coponent(1) --> line 1 in the tensor. 
	// - iUEnergy.coponent(2) --> line 2 in the tensor. 
	tmp<volTensorField> calculateEnergyFlux(tmp<volTensorField> ioFluxTensorPtr, volVectorField& iUEnergy);

	// - Multiplies the tensor flux by iUEnergy. 
	// - iUEnergy.coponent(0) --> line 0 in the tensor. 
	// - iUEnergy.coponent(1) --> line 1 in the tensor. 
	// - iUEnergy.coponent(2) --> line 2 in the tensor. 
	tmp<volTensorField> calculateEnergyFlux(surfaceScalarField& iphi, volVectorField& iU, volVectorField& iUEnergy);

	// - Return a new field with the value of Uc in every cell. 
	void setUc(volVectorField& iU); 



	// ==========================================================
	// ========================================================== util. 
	// ==========================================================

	scalar Integrate(const volScalarField & field) { 
			scalar sfield;
			sfield = sum(field*zoneSelector).value(); 
			reduce(sfield,sumOp<scalar>());
			return sfield/ZoneVolume;
	}

	scalar Integrate(const tmp<volScalarField> & field) { 
			scalar sfield;
			sfield = sum(field*zoneSelector).value(); 
			reduce(sfield,sumOp<scalar>());
			return sfield/ZoneVolume;
	}

	vector Integrate(const volVectorField & field) { 
			vector sfield;
			sfield = sum(field*zoneSelector).value(); 
			reduce(sfield,sumOp<vector>());
			return sfield/ZoneVolume;
	}


	vector Integrate(const tmp<volVectorField> & field) { 
			vector sfield;
			sfield = sum(field*zoneSelector).value(); 
			reduce(sfield,sumOp<vector>());
			return sfield/ZoneVolume;
	}

	tensor Integrate(const tmp<volTensorField> & field) { 
			tensor sfield;
			sfield = sum(field*zoneSelector).value(); 
			reduce(sfield,sumOp<tensor>());
			return sfield/ZoneVolume;
	}

	tensor Integrate(const volTensorField& field) { 
			tensor sfield;
			sfield = sum(field*zoneSelector).value(); 
			reduce(sfield,sumOp<tensor>());
			return sfield/ZoneVolume;
	}

	// =======================================================================

	// Checking the momentum balance of this timestep. 
	void checkMomentumBalance_Timestep();
	

};
#endif
