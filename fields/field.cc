// gemc headers
#include "field.h"
#include "string_utilities.h"
#include "fieldFactory.h"
#include "multipoleField.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// G4 headers
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4CashKarpRKF45.hh"
#include "G4SimpleHeum.hh"
#include "G4SimpleRunge.hh"
#include "G4NystromRK4.hh"
#include "G4CachedMagneticField.hh"

// this class serves as a dispatcher for the
// various format of magnetic fields
// availiable formats:
// - simple (uniform fields)
// - map
void gfield::create_MFM()
{

	// fields can be uniform, mapped
	if(format == "simple" && symmetry == "uniform")
		create_simple_MFM();

	// fields can be multipole, mapped
	if(format == "simple" && symmetry == "multipole")
		create_simple_multipole_MFM();

	if(format == "map")
	{
		fFactory->loadFieldMap(map, verbosity);

		G4Mag_UsualEqRhs*       iEquation    = new G4Mag_UsualEqRhs(map);
		G4MagIntegratorStepper* iStepper     = createStepper(integration, 	iEquation);
		G4ChordFinder*          iChordFinder = new G4ChordFinder(map, minStep, iStepper);

		// caching does not seem to help for dipole-y
		// will it help for other field maps?
		G4MagneticField *pCachedMagField = new G4CachedMagneticField(map, 1 * m);
		MFM = new G4FieldManager(pCachedMagField, iChordFinder);

		G4double minEps = 0.1;  //   Minimum & value for smallest steps
		G4double maxEps = 1.0;  //   Maximum & value for largest steps
		
		MFM->SetMinimumEpsilonStep( minEps );
		MFM->SetMaximumEpsilonStep( maxEps );
 		MFM->SetDeltaOneStep(0.01 * mm);
		MFM->SetDeltaIntersection(0.01 * mm);
	}

}

void gfield::create_simple_MFM()
{
	vector < string > dim = getStringVectorFromString(dimensions);
	
	if (dim.size() != 3)
		cout << "   !!! Error: dimension of field " << name << " are wrong: "
	  		 << dimensions << " has dimension " << dim.size() << endl;

	const G4ThreeVector constField(get_number(dim[0]), get_number(dim[1]), get_number(dim[2]));
	G4UniformMagField*      magField     = new G4UniformMagField(constField);

	G4Mag_UsualEqRhs*       iEquation    = new G4Mag_UsualEqRhs(magField);
	G4MagIntegratorStepper* iStepper     = createStepper(integration, iEquation);
	G4ChordFinder*          iChordFinder = new G4ChordFinder(magField, minStep,	iStepper);

	MFM = new G4FieldManager(magField, iChordFinder);

	if (verbosity > 1)
	{
		cout << "  >  <" << name << ">: uniform magnetic field is built." << endl;
	}
}

void gfield::create_simple_multipole_MFM()
{
	vector < string > dim = getStringVectorFromString(dimensions);
	
	if (dim.size() != 7)
		cout << "   !!! Error: dimension of field " << name << " are wrong: "
				<< dimensions << " has dimension " << dim.size() << endl;

	multipoleField* magField = new multipoleField(atoi(dim[0].c_str()), get_number(dim[1]), get_number(dim[2]), get_number(dim[3]),
			                                     get_number(dim[4]), get_number(dim[5]), dim[6]);
	
	
	G4Mag_UsualEqRhs* iEquation      = new G4Mag_UsualEqRhs(magField);
	G4MagIntegratorStepper* iStepper = createStepper(integration, iEquation);
	G4ChordFinder* iChordFinder      = new G4ChordFinder(magField, minStep, iStepper);

	MFM = new G4FieldManager(magField, iChordFinder);

	if (verbosity > 1)
	{
		cout << "  >  <" << name << ">: multipole magnetic field is built with "
			 << dimensions << endl;
	}
}

G4MagIntegratorStepper *createStepper(string sname, G4Mag_UsualEqRhs* ie)
{
	if (sname == "G4CashKarpRKF45")	      return new G4CashKarpRKF45(ie);
	if (sname == "G4ClassicalRK4")	      return new G4ClassicalRK4(ie);
	if (sname == "G4SimpleHeum")		      return new G4SimpleHeum(ie);
	if (sname == "G4SimpleRunge")	         return new G4SimpleRunge(ie);
	if (sname == "G4ImplicitEuler")	      return new G4ImplicitEuler(ie);
	if (sname == "G4ExplicitEuler")	      return new G4ExplicitEuler(ie);
	if (sname == "G4HelixImplicitEuler")	return new G4HelixImplicitEuler(ie);
	if (sname == "G4HelixExplicitEuler")	return new G4HelixExplicitEuler(ie);
	if (sname == "G4HelixSimpleRunge")	   return new G4HelixSimpleRunge(ie);
	if (sname == "G4NystromRK4")	         return new G4NystromRK4(ie);

	// if requested is not found return NULL
	cout << "  !!! Error: stepper " << sname << " is not defined " << endl;
	return NULL;
}

void gfield::initialize(goptions Opt)
{
	string hd_msg = "  >> fields Init: ";

	vector<aopt> FIELD_SCALES_OPTION = Opt.getArgs("SCALE_FIELD");
	for (unsigned int f = 0; f < FIELD_SCALES_OPTION.size(); f++)
	{
		vector < string > scales = getStringVectorFromStringWithDelimiter(FIELD_SCALES_OPTION[f].args, ",");
		if(scales.size() == 2)
		{
			if (scales[0].find(name) != string::npos)
				scaleFactor = get_number(scales[1]);
		}
	}

	vector<aopt> FIELD_PROPERTIES = Opt.getArgs("FIELD_PROPERTIES");
	for (unsigned int f = 0; f < FIELD_PROPERTIES.size(); f++)
	{
		vector < string > attributes = getStringVectorFromStringWithDelimiter(FIELD_PROPERTIES[f].args, ",");
		if(attributes.size() > 2)
		{
			if(attributes[0].find(name) != string::npos)
			{
				minStep = get_number(attributes[1]);
				integration = trimSpacesFromString(attributes[2]);

				if(map)
					if(attributes.size() == 4)
						map->interpolation = trimSpacesFromString(attributes[3]);
			}
		}
	}


	if(map)
	{
		map->scaleFactor = scaleFactor;
		map->initializeMap();
		map->verbosity = verbosity;
	}
}

///< Overloaded "<<" for gfield class. Dumps infos on screen.
ostream &operator<<(ostream &stream, gfield gf)
{
	cout << "  > Field name:           "   << gf.name << endl;
	cout << "    - factory:            "   << gf.factory << endl;
	cout << "    - comment:            "   << gf.description << endl;
	cout << "    - format:             "   << gf.format << endl;
	cout << "    - symmetry:           "   << gf.symmetry << endl;
	cout << "    - scale factor:       "   << gf.scaleFactor << endl;
	cout << "    - integration method: "   << gf.integration << endl;
	cout << "    - minimum Step:       "   << gf.minStep << " mm" << endl;
	if (gf.dimensions != "na" && gf.format == "simple")
		cout << "    - dimensions:         " << gf.dimensions << endl;

	if (gf.dimensions == "na" && gf.format == "map")
	{
		cout << "    - map identifier:     " << gf.map->identifier << endl;

		for (unsigned int i = 0; i < gf.map->coordinates.size(); i++)
			cout << "    - Coordinate:         " << gf.map->coordinates[i];
		
		cout << "    - Map Field Unit:     " << gf.map->unit << endl;
		cout << "    - Map Interpolation:  " << gf.map->interpolation << endl;

		cout << "    - map origin:         x=" << gf.map->mapOrigin[0]
			 << "mm, y=" << gf.map->mapOrigin[1]
		     << "mm, z=" << gf.map->mapOrigin[2] << "mm" << endl;
	}
	return stream;
}




