
	globalMeshData globalmeshdata(mesh); 
	const globalIndex& globalindexer = globalmeshdata.globalPointNumbering();

	Info << " Local " << globalindexer.isLocal(401) << endl;//<< " " << globalindexer.toLocal(400) << endl;


	// Building the zone selection idicator field. 
	volScalarField zoneSelector
	(
            IOobject
            (
                    "ZoneSelector",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
	);

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
	zoneSelector.write();



	scalar ZoneVolume = sum(zoneSelector).value(); 
	Pout << "CPU: " << ZoneVolume << endl;
	reduce(ZoneVolume,sumOp<scalar>());

	Info << "The zone volume is " << ZoneVolume << endl; 

