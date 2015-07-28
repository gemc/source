/// \file field.h 
/// Defines the gemc Field class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef FIELD_H
#define FIELD_H 1



// The magnetic fields are self-descriptive, using an XML reader
// The maps are stored in several format. 
// Each format corresponds to a gfield factory. Example: ASCII, EVIO.
//
// At start time gemc loads a map<string, string> looking for 
// eligible field definitions. Value is the factory type.
// Each gfield factory will then produce a map<string, gfield> 
// from this list, based on the file format.
//
// The gfield is loaded in the map if the header of 
// the file is understood and provides a valid gfield description
//
// At detector construction time the mfield entry in detector.h, if 
// different than "no" must be associated to a key of the  map<string, gfield>.
//
// At that point the MFM is checked on the gfield. If NULL, the MFM is created. Either way, the MFM will be 
// associated with the logical volume.

// gemc headers
#include "options.h"
#include "mappedField.h"
// forward declaration of fieldFactory
class fieldFactory;

// C++ headers
#include <string>
using namespace std;

// G4 headers
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


/// \class gfield
/// <b>gfield </b>\n\n
/// Defines the gemc Field class.\n
/// The class parameters are filled from the field factory
class gfield
{
	public:
		gfield(){;}
		gfield(goptions opts)
		{
			// initialize Magnetic Field Manager and Mapped field to NULL
			MFM         = NULL;
			map         = NULL;
			symmetry    = "na";
			format      = "na";
			dimensions  = "na";
			scaleFactor	= 1;
			minStep     = 0.1*mm;
			integration = "ClassicalRK4";
			verbosity   = opts.optMap["FIELD_VERBOSITY"].arg;
		}
	 ~gfield(){}
	
	public:
		string name;            ///< Field name - used as key in the map<string, gfield>
		string symmetry;        ///< Field symmetry
		string format;          ///< Field format (available: simple (for uniform) and map)
		string factory;         ///< Field factory (format of magnetic field)
		string description;     ///< Field Description
		string dimensions;      ///< Field dimensions (with units), for non-mapped fields
		string integration;     ///< Integration Method
		double verbosity;       ///< Log verbosity
		double minStep;         ///< Minimum Step for the G4ChordFinder
		string unit;            ///< Field Unit
	
		// Scale factor is set from options
		double scaleFactor;     
		void initialize(goptions);
	
		// creates simple magnetic field manager (uniform fields, etc)
		void create_simple_MFM();
		void create_simple_multipole_MFM();
	
		// mapped Field. We need to factory to load the map
		gMappedField *map;       ///< Mapped Field
		fieldFactory *fFactory;  ///< fieldFactory that created the field
	
	private:
		G4FieldManager *MFM;             	///< G4 Magnetic Field Manager
		void create_MFM();                ///< Creates the G4 Magnetic Field Manager
	
	public:
		// Returns Magnetic Field Manager Pointer
		// creates one if it doesn't exist
		G4FieldManager* get_MFM()
		{
			if(MFM == NULL)
				create_MFM();
				
			return MFM;
		} 	

		///< Overloaded "<<" for gfield class. Dumps infos on screen.
		friend ostream &operator<<(ostream &stream, gfield gf);
				
};


// creates custom integrator based on the equation
G4MagIntegratorStepper *createStepper(string sname, G4Mag_UsualEqRhs *iE); 


#endif
