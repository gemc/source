#ifndef BONUS_HITPROCESS_H
#define BONUS_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition
/// \class bonus_HitProcess
/// <b> Ceontral Time of Flight Hit Process Routine</b>\n\n
/// The Calibration Constants are:\n
/// - VEF is the effective velocity of propogation in the scintillator

class bonus_HitProcess : public HitProcess
{
	public:
	
		~bonus_HitProcess(){;}
	
		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);
				
		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
		// creates the HitProcess
		static HitProcess *createHitClass() {return new bonus_HitProcess;}
};

#endif
