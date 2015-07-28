#ifndef multipoleField_HH
#define multipoleField_HH

// gemc headers
#include "string_utilities.h"

// G4 headers
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

class multipoleField: public G4MagneticField
{
	public:

		multipoleField(int Npole, G4double scale,
					   G4double x, G4double y, G4double z,  // origin
 				   	G4double rot, string ROTaxis);       // axis of rotation

		virtual ~multipoleField();
	
		virtual void GetFieldValue(const G4double pos[4], G4double *MagField) const;

		int polenumber;
		G4double origin[3];
		G4double strength;
		G4double rotation;
		string rotaxis;

};

#endif
