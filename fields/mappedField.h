/// \file field.h 
/// Defines the gemc Field class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef MAPPED_FIELD_H
#define MAPPED_FIELD_H 1

// G4 headers
#include "G4MagneticField.hh"


// c++ headers
#include <vector>
#include <string>
#include <iostream>

using namespace std;


// defines one dimension of the field
// number of points, range and units
class gcoord
{
public:
	gcoord(string nm, unsigned int n, double m, double M, string u, int s)
	{
		name  = nm;
		np    = n;
		min   = m;
		max   = M;
		unit  = u;
		speed = s;
	}
	
public:
	string       name;
	unsigned int np;
	double       min;
	double       max;
	string       unit;
	int          speed;     // 0 is the slowest varying coordinate
	
	friend ostream &operator<<(ostream &stream, gcoord gc);
};



// define a mapped field
/// \class gMappedField
/// <b>gMappedField </b>\n\n
/// This class defines gemc Mapped Electro-Magnetic Fields.\n
/// It implements G4MagneticField::GetFieldValue
/// that returns a magnetic field value at a point in space
class gMappedField : public G4MagneticField
{
public:
	gMappedField(string i, string s): identifier(i), symmetry(s) {
		mapOrigin[0]  = 0;
		mapOrigin[1]  = 0;
		mapOrigin[2]  = 0;
		mapRotation[0] = 0;
		mapRotation[1] = 0;
		mapRotation[2] = 0;
		scaleFactor   = 1;
		unit          = "gauss";
		interpolation = "linear";
		verbosity     = 0;
	}
	~gMappedField(){;}
	
	string identifier;          ///< Pointer to map in factory (for example, hostname / filename with path / date)
	string symmetry;            ///< map symmetry
	vector<gcoord> coordinates; ///< Vector size depend on the symmetry

	// set by gfield::initialize
	double mapOrigin[3];        ///< Displacement of map. This is used in GetFieldValue
	double mapRotation[3];      ///< Rotation of map. This is used in GetFieldValue
	double scaleFactor;         ///< copy of the gfield scaleFactor

	string unit;                ///< field unit in the map
	string interpolation;       ///< map interpolation technique. Choices are "none", "linear", "quadratic"
	int verbosity;              ///< map verbosity
	
	// field depending on 3D map
	double ***B1_3D;
	double ***B2_3D;
	double ***B3_3D;
	
	// field depending on 2D map
	double **B1_2D;
	double **B2_2D;
	
	// these are initialized based on the map
	// symmetry and coordinates
	// it will avoid the time spend in retrieving
	// the infos in GetFieldValue
	double *startMap;
	double *endMap;	
	double *cellSize;
	unsigned int *np;
	void initializeMap();
	
	gcoord getCoordinateWithSpeed(int speed);   ///< return coordinate based on speed
	gcoord getCoordinateWithName(string name);  ///< return coordinate based on type
	
	// returns the field at point x. This is a dispatcher for the various symmetries below
	void GetFieldValue( const double x[3], double *Bfield) const;
	
	void GetFieldValue_Dipole( const double x[3], double *Bfield, int FIRST_ONLY) const;
	void GetFieldValue_Cylindrical( const double x[3], double *Bfield, int FIRST_ONLY) const;
	void GetFieldValue_phiSegmented( const double x[3], double *Bfield, int FIRST_ONLY) const;
	void GetFieldValue_cartesian3d( const double x[3], double *Bfield, int FIRST_ONLY) const;
	
	// precalculating values of the rotation angles so we don't do it at GetFieldValue time
	double sinAlpha, cosAlhpa;
	double sinBeta, cosBeta;
	double sinGamma, cosGamma;
	void RotateField( double *Bfield) const;

	// we want to rotate the field (axes), not the point
	// so each rotation is the inverse of the point rotation
	inline const double yRotX(double p[3]) const {return  p[1]*cosAlhpa + p[2]*sinAlpha;}
	inline const double zRotX(double p[3]) const {return -p[1]*sinAlpha + p[2]*cosAlhpa;}

	inline const double xRotY(double p[3]) const {return p[0]*cosBeta - p[2]*sinBeta;}
	inline const double zRotY(double p[3]) const {return p[0]*sinBeta + p[2]*cosBeta;}

	inline const double xRotZ(double p[3]) const {return  p[0]*cosGamma + p[1]*sinGamma;}
	inline const double yRotZ(double p[3]) const {return -p[0]*sinGamma + p[1]*cosGamma;}

};

#endif







