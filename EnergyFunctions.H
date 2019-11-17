/*
scalar IntegrateVector(const volScalarField & field) { 
		scalar sfield = sum(field).value(); 
		reduce(sfield,sumOp<scalar>());
		return sfield;
}

scalar IntegrateVector(const tmp<volScalarField> & field) { 
		scalar sfield = sum(field).value(); 
		reduce(sfield,sumOp<scalar>());
		return sfield;
}

vector IntegrateVector(const volVectorField & field) { 
		vector sfield = sum(field).value(); 
		reduce(sfield,sumOp<vector>());
		return sfield;
}


vector IntegrateVector(const tmp<volVectorField> & field) { 
		vector sfield = sum(field).value(); 
		reduce(sfield,sumOp<vector>());
		return sfield;
}

vector IntegrateVector(const tmp<volTensorField> & field) { 
		tensor sfield = sum(field).value(); 
		reduce(sfield,sumOp<tensor>());
		return vector(sfield.xx(),sfield.yy(),sfield.zz());
}

tensor IntegrateTensor(const tmp<volTensorField> & field) { 
		tensor sfield = sum(field).value(); 
		reduce(sfield,sumOp<tensor>());
		return sfield;
}

tensor IntegrateTensor(const volTensorField& field) { 
		tensor sfield = sum(field).value(); 
		reduce(sfield,sumOp<tensor>());
		return sfield;
}

*/

