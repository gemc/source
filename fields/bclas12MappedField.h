/// \file field.h 
/// Defines the gemc Field class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef CLAS12BINARYMAPPED_FIELD_H
#define CLAS12BINARYMAPPED_FIELD_H 1

// G4 headers
#include "G4MagneticField.hh"

extern "C" {
#include "magfield.h"
#include "magfieldio.h"
#include "munittest.h"
#include "magfieldutil.h"
}

// c++ headers
#include <vector>
#include <string>
#include <iostream>
using namespace std;



// define a mapped field
/// \class gclas12BinaryMappedField
/// <b>gclas12BinaryMappedField </b>\n\n
/// This class defines gemc gclas12BinaryMappedField .\n
/// It implements G4MagneticField::GetFieldValue
/// that returns a magnetic field value at a point in space

class gclas12BinaryMappedField : public G4MagneticField
{
public:
	gclas12BinaryMappedField(string identity): identifier(identity) {

		// initialize them to zero here, can be set by derived gfield::initialize
		mapOrigin[0]   = 0;
		mapOrigin[1]   = 0;
		mapOrigin[2]   = 0;
		mapRotation[0] = 0;
		mapRotation[1] = 0;
		mapRotation[2] = 0;
		scaleFactor    = 1;
		interpolation = "linear";
		verbosity      = 0;

	}
	~gclas12BinaryMappedField(){;}
	
	int verbosity;              ///< map verbosity
	string identifier;          ///< Pointer to map in factory (for example, hostname / filename with path / date)

	// set by gfield::initialize
	double mapOrigin[3];        ///< Displacement of map. This is used in GetFieldValue
	double mapRotation[3];      ///< Rotation of map. This is used in GetFieldValue
	double scaleFactor;         ///< copy of the gfield scaleFactor
	string interpolation;       ///< map interpolation technique. Choices are "none", "linear", "quadratic"
	string unit;                ///< field unit in the map

	// returns the field at point x. This is a dispatcher for the various symmetries below
	void GetFieldValue( const double x[3], double *Bfield) const;

	// fields filenames
	string solenoidMapFileName;
	string torusMapFileName;

	MagneticFieldPtr solenoidPtr;
	MagneticFieldPtr torusPtr;

	// map pointers
	// Uses David's routine to return the BX BY BZ components
	FieldValuePtr combinedValuePtr;

	// precalculating values of the rotation angles so we don't do it at GetFieldValue time
	double sinAlpha, cosAlhpa;
	double sinBeta, cosBeta;
	double sinGamma, cosGamma;

	// sets the values above
	void initializeMap();

	void RotateField( double *Bfield) const;

	// we want to rotate the field (axes), not the point
	// so each rotation is the inverse of the point rotation
	inline const double yRotX(double p[3]) const {return  p[1]*cosAlhpa + p[2]*sinAlpha;}
	inline const double zRotX(double p[3]) const {return -p[1]*sinAlpha + p[2]*cosAlhpa;}

	inline const double xRotY(double p[3]) const {return p[0]*cosBeta - p[2]*sinBeta;}
	inline const double zRotY(double p[3]) const {return p[0]*sinBeta + p[2]*cosBeta;}

	inline const double xRotZ(double p[3]) const {return  p[0]*cosGamma + p[1]*sinGamma;}
	inline const double yRotZ(double p[3]) const {return -p[0]*sinGamma + p[1]*cosGamma;}

	
	void defineNamesAndType(string field_dir);

private:


};

#endif














