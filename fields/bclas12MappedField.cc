// gemc include
#include "bclas12MappedField.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


void gclas12BinaryMappedField::GetFieldValue(const double x[3], double *bField) const
{
	static int FIRST_ONLY;
	
	// displacement point
	double rpoint[3] = {x[0] - mapOrigin[0], x[1] - mapOrigin[1], x[2] - mapOrigin[2]};
	

	bField[0] = bField[1] = bField[2] = 0;


	// Uses David's routine to return the BX BY BZ components
	FieldValuePtr combinedValuePtr = (FieldValuePtr) malloc(sizeof (FieldValue));
	//getCompositeFieldValue(combinedValuePtr, x[0], x[1], x[2], torus, solenoid);//torus and solenoid were declared in the header file

	bField[0] = combinedValuePtr->b1;
	bField[1] = combinedValuePtr->b2;
	bField[2] = combinedValuePtr->b3;

	RotateField(bField);

	// we don't worry about computer speed
	// if verbosity is set this high
	// so we can output units as well
	// add any useful information here
	if(verbosity>3 && FIRST_ONLY != 99) {
		cout << "  > Track position in magnetic field map, with displacement and rotations (x,y,z)/cm:"
		<< "("  << x[0]/cm << ", "
		<< x[1]/cm << ", "
		<< x[2]/cm << ") cm,  " << endl;
		cout << "    Cylindrical: ";
		cout << "loc. pos. = ("    << x[0]/cm << ", " << x[1]/cm << ", " << x[2]/cm << ") cm,  ";
		cout << "B = ("   << bField[0]/gauss << ",  " << bField[1]/gauss << ",  " << bField[2]/gauss << ") gauss " << endl;
	}

	if(verbosity == 99) FIRST_ONLY = 99;
	
}

void gclas12BinaryMappedField::RotateField( double *Bfield) const  {

	// rotating the fields
	if(mapRotation[0] != 0) {
		double yPrime = yRotX(Bfield);
		double zPrime = zRotX(Bfield);
		Bfield[1] = yPrime;
		Bfield[2] = zPrime;
	}

	if(mapRotation[1] != 0) {
		double xPrime = xRotY(Bfield);
		double zPrime = zRotY(Bfield);
		Bfield[0] = xPrime;
		Bfield[2] = zPrime;
	}

	if(mapRotation[2] != 0) {
		double xPrime = xRotZ(Bfield);
		double yPrime = yRotZ(Bfield);
		Bfield[0] = xPrime;
		Bfield[1] = yPrime;
	}
}




void gclas12BinaryMappedField::initializeMap()
{
	// setting rotation sin and cosines
	sinAlpha = sin(mapRotation[0]);
	cosAlhpa = cos(mapRotation[0]);
	sinBeta = sin(mapRotation[1]);
	cosBeta = cos(mapRotation[1]);
	sinGamma = sin(mapRotation[2]);
	cosGamma = cos(mapRotation[2]);
}









