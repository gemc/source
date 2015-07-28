#ifndef MDetectorMessenger_h
#define MDetectorMessenger_h 1

// G4 headers
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

// gemc headers
#include "MDetectorConstruction.h"

// Class definition
class MDetectorMessenger : public G4UImessenger
{
	public:
		MDetectorMessenger(MDetectorConstruction*);
		~MDetectorMessenger();
		void SetNewValue(G4UIcommand*, G4String);
		
	public:
		
	private:
		MDetectorConstruction   *MDetector;
		G4UIdirectory           *gemcDir;
		G4UIcmdWithoutParameter *UpdateGeoCmd;
		
};

#endif
