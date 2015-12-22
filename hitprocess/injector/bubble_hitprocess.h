#ifndef bubble_HITPROCESS_H
#define bubble_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition
/// \class bubble_HitProcess
/// <b> Bubble Hit Process Routine</b>\n\n
/// The Calibration Constants are:\n
/// - VEF is the effective velocity of propogation in the scintillator

class bubble_HitProcess : public HitProcess
{
	public:
	
		~bubble_HitProcess(){;}
	
		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);
				
		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
		// creates the HitProcess
		static HitProcess *createHitClass() {return new bubble_HitProcess;}
};

#endif
