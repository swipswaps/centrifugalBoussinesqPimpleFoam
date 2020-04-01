#ifndef EnergyCalculations_H
#define EnergyCalculations_H 

#include "fvCFD.H"
#include "timeSelector.H"
#include "fvc.H"
#include "fvCFD.H"

#include "meshSearch.H"


//#define TESTS


#define def_fields(type,term,ctype)\
	vol##type##Field mean_##ctype##_##term; \
	vol##type##Field perb_##ctype##_##term;\
	vol##type##Field total_##ctype##_##term;

	



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
	def_fields(Scalar,dUdt,energy);    		// dt terms. 
	def_fields(Tensor,diffusion,energy);    	// diffusion terms.
	def_fields(Tensor,gradUsqr,energy);    	// The diffusion terms after the integration. 
	def_fields(Scalar,pressure,energy);    	// pressure terms.
	def_fields(Scalar,nudging,energy);    	// nudging terms.
	def_fields(Scalar,potential,energy);    	// potential energy terms.        

	// - Convection terms. 
	volTensorField energy_fullFlux; 	// total flux of energy. 
	volTensorField mean_meanMeanFlux; 	// mean flux of mean kinetic energy. 
	volTensorField mean_perbPerbFlux;       // = Conversion term + reynolds flux. 
	volTensorField perb_perbMeanFlux;       // = Conversion term (Method II). 
	volTensorField perb_meanPerbFlux; 	// = Mean Flux of perturbation kinetic energy. 
	volTensorField perb_perbPerbFlux; 	// = Perturbation flux of perturbation kinetic energy. 

	volTensorField mean_perb_conversion; 	// = Direct calculation of the conversion term. (Method I).


	def_fields(Tensor,diffusion,momentum);    	// diffusion terms.
	def_fields(Vector,diffusion,eqn);		// the diffusion terms in each wind component. 
	def_fields(Vector,diffusion,eqn_energy);        // the energy diffusion terms in each wind component. 

#ifdef TESTS
	volScalarField energy_fullFlux_eqn;   // U&fvc::div(phi,U);
	volScalarField mean_meanMeanFlux_eqn; // meanU&fvc::div(meanphi,meanU);
	volScalarField mean_perbPerbFlux_eqn; // tagU&fvc::div(tag_phi,meanU);
	volScalarField perb_meanPerbFlux_eqn; // tagU&fvc::div(meanphi,tagU);   
	volScalarField perb_perbMeanFlux_eqn; // tagU&fvc::div(tag_phi,meanU);
	volScalarField perb_perbPerbFlux_eqn; // tagU&fvc::div(tag_phi,tagU);    


	// should be 0. 
	volScalarField perb_meanMeanFlux_eqn;
	volScalarField mean_meanPerbFlux_eqn;  
	volScalarField mean_perbMeanFlux_eqn;  

	// DEBUGGING FIELDS. 
	volScalarField mean_fullFlux_eqn;
	volScalarField perb_fullFlux_eqn;
	volScalarField mean_fullTagFlux_eqn;
	volScalarField mean_fullMeanFlux_eqn;
	
	volScalarField mean_meanFullFlux_eqn; 
	volScalarField mean_perbFullFlux_eqn; 
	volVectorField meanPerbFlux_eqn;
#endif
	



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

	// Calculates the correct term-wise laplacian of the field U. 
	tmp<volTensorField> calculategradKgrad(volVectorField& U);


	// ==========================================================
	// ========================================================== tests. 
	// ==========================================================

	// test the mean and perturbation terms; 
	void test_energy_dUdt(); 
	void test_energy_Pressure(); 
	void test_energy_Diffusion();
	void test_energy_Convection();
	void test_energy_Nudging(); 
	void test_enetgy_Potential();


#ifdef TESTS
	void update_energy_Convection_eqn_sum();

	void test_energy_Convection_Mean_perbperb(tmp<volTensorField> mean_perbPerbFlux_current);
	void test_energy_Convection_perb_perbMeanFlux(tmp<volTensorField> perb_perbMeanFlux_current);
	void test_energy_Convection_perb_meanPerbFlux(tmp<volTensorField> perb_meanPerbFlux_current);
	void test_energy_Convection_perb_perbPerbFlux(tmp<volTensorField> perb_perbPerbFlux_current);
	void test_energy_Convection_eqn();
#endif


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
