/// \file MSteppingAction.h
/// Defines the gemc Stepping Action class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef MSteppingAction_h
#define MSteppingAction_h 1

// G4 headers
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

// gemc headers
#include "options.h"

class MSteppingAction : public G4UserSteppingAction
{
	public:
		MSteppingAction(goptions);
		virtual ~MSteppingAction();
		
		goptions gemcOpt;            ///< gemc options map
		double energyCut;            ///< Set to the ENERGY_CUT from options. This avoids a billion lookups in the map.
		double max_x_pos;            ///< Max X Position in millimeters.
		double max_y_pos;            ///< Max Y Position in millimeters.
		double max_z_pos;            ///< Max Z Position in millimeters.
				
		// checking if track get stuck.
		// if after 10 times the oldpos is the same as new pos,
		// the track will be killed
//		G4ThreeVector oldpos;
//		int            nsame;
		
		void UserSteppingAction(const G4Step*);
};

#endif
