    Info<< "Reading thermophysical properties\n" << endl;

    Info<< "Reading field T mean\n" << endl;
    volScalarField Tmean
    (
        IOobject
        (
            "Tmean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T (actually T perturbation)\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField CentrifugalFactor
    (
        IOobject
        (
            "CentrifugalFactor",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    dimensionedScalar TRef = sum(Tmean*mesh.V())/sum(mesh.V());

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U mean\n" << endl;
    volVectorField Umean
    (
        IOobject
        (
            "Umean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U (actuallty U perturbation)\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );




    Info<< "Creating field Utotal = U (perturbation) + Umean\n" << endl;
    volVectorField Utotal
    (
        IOobject
        (
            "Utotal",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U+Umean
    );
    Utotal.write(); 

    Info<< "Creating field Ttotal = T (perturbation) + Tmean\n" << endl;
    volScalarField Ttotal
    (
        IOobject
        (
            "Ttotal",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T+Tmean
    );
    Ttotal.write(); 

    volVectorField Uwhitenoise
    (
    	IOobject
        (
            "Uwhitenoise",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
     );

    volScalarField Twhitenoise
    (
    	IOobject
        (
            "Twhitenoise",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh
     );
     
   
     
     
     

    Info<< "Reading/calculating face flux field of mean U: phimean\n" << endl;

    surfaceScalarField phimean
    (
	    IOobject
	    (
		"phimean",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	    ),
	    linearInterpolate(Umean) & mesh.Sf()
    );


    #include "createPhi.H"

    #include "readTransportProperties.H"

    Info << "Calculating the mean rho" << nl; 

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );


    // Kinematic density for buoyancy force
    volScalarField rhok_tag
    (
        IOobject
        (
            "rhok_tag",
            runTime.timeName(),
            mesh
        ),
        - beta*T
    );
    
    // kinematic turbulent thermal thermal conductivity m2/s
    Info<< "Reading field alphat\n" << endl;
    volScalarField alphat
    (
        IOobject
        (
            "alphat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

    volTensorField AnisotropicDiffusion
    (
	IOobject
    	(
	"AnisotropicDiffusion",
	runTime.timeName(),
	mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
	),
	mesh
    );

