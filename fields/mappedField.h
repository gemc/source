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

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// defines one dimension of the field
// number of points, range and units
class gcoord
{
	public:
		gcoord(string nm = "na", unsigned int n = 0, double m = 0, double M=0, string u="na", int s=0) :
			name  ( nm ), np    ( n  ), min   ( m  ), max   ( M  ), unit  ( u  ), speed ( s  )
		{ }
	
	public:
		string       name;
		unsigned int np;
		double       min;
		double       max;
		string       unit;
		int          speed;     // 0 is the slowest varying coordinate



		friend ostream &operator<<(ostream &stream, gcoord gc)
		{
			cout << gc.name  << ": np="    << gc.np ;
			
			if(gc.unit == "mm" || gc.unit == "m" ||  gc.unit == "cm")
			{
				cout << ", min="   << gc.min/cm << " cm"
				     << ", max="   << gc.max/cm << " cm";
			}
			if(gc.unit == "deg" || gc.unit == "rad" )
			{
				cout << ", min="   << gc.min/degree << " deg"
				     << ", max="   << gc.max/degree << " deg";
			}
			
			cout << ", index speed=" << gc.speed << endl;
			return stream;
		}
};


// define a mapped field
/// \class gMappedField
/// <b>gMappedField </b>\n\n
/// This class defines gemc Mapped Electro-Magnetic Fields.\n
/// The function G4MagneticField function GetFieldValue
/// returns a magnetic field value at a point in space
class gMappedField : public G4MagneticField
{
   public:

      enum class MapInterpolation {none, linear, quadratic}; 
      enum class SymmetryAxis {x, y, z}; 

      gMappedField(string i, string s)
      {		
         symmetryAxis  = SymmetryAxis::z;
         symmetry      = s;
         identifier    = i;
         mapOrigin[0]  = 0;
         mapOrigin[1]  = 0;
         mapOrigin[2]  = 0;
         scaleFactor   = 1;
         unit          = "gauss";
         interpolation = MapInterpolation::none;
         verbosity     = 0;
         func_ptr_GetFieldValue = 0;

         double sqrt3 = sqrt(3.0);
         //     angle=     0         60        120   180         240         300   360 
         segment_sin = { 0.0, sqrt3/2.0, sqrt3/2.0,  0.0, -sqrt3/2.0, -sqrt3/2.0,  0.0 };
         segment_cos = { 1.0,       0.5,      -0.5, -1.0,       -0.5,        0.5,  1.0 };

      }
      ~gMappedField(){;}

      SymmetryAxis  symmetryAxis;
      string symmetry;            ///< map symmetry
      string identifier;          ///< Pointer to map in factory (for example, hostname / filename with path / date)
      vector<gcoord> coordinates; ///< Vector size depend on the symmetry

      double mapOrigin[3];        ///< Displacement of map. This is used in GetFieldValue
      double scaleFactor;         ///< copy of the gfield scaleFactor
      string unit;                ///< field unit in the map
      MapInterpolation interpolation;
      //string interpolation;       ///< map interpolation technique. Choices are "none", "linear", "quadratic"
      int verbosity;              ///< map verbosity

      vector<double>   segment_sin;
      vector<double>   segment_cos;

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
      double *cellSize;
      unsigned int *np;
      void initializeMap();

      const gcoord& getCoordinateWithSpeed(int speed);   ///< return coordinate based on speed
      gcoord getCoordinateWithName(string name);  ///< return coordinate based on type
      gcoord dummy;

      // returns the field at point x. This is a dispatcher for the various symmetries below
      void GetFieldValue( const double x[3], double *Bfield) const;

      void SetFunctionPointer() const;
      mutable void (gMappedField::*func_ptr_GetFieldValue)( const double x[3], double *Bfield, int) const ;

      void GetFieldValue_Dipole( const double x[3], double *Bfield, int FIRST_ONLY) const;
      void GetFieldValue_Cylindrical( const double x[3], double *Bfield, int FIRST_ONLY) const;
      void GetFieldValue_phiSegmented( const double x[3], double *Bfield, int FIRST_ONLY) const;

};

#endif







