#ifndef DC_HITPROCESS_H
#define DC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// Class definition
class dc_HitProcess : public HitProcess
{
	public:
	
		~dc_HitProcess(){;}
	
		// initialize private vars
		 dc_HitProcess()
		{
			mini_stagger_shift_R2 = gemcOpt.optMap["DC_MSTAG_R2"].arg*mm;   // Mini Stagger shift for Region 2
			mini_stagger_shift_R3 = gemcOpt.optMap["DC_MSTAG_R3"].arg*mm;   // Mini Stagger shift for Region 3
		
			NWIRES = 113;
		}
	
		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);
		
		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
		// creates the HitProcess
		static HitProcess *createHitClass() {return new dc_HitProcess;}


	private:
		double mini_stagger_shift_R2;   // Mini Stagger shift for Region 2
		double mini_stagger_shift_R3;   // Mini Stagger shift for Region 3

		double NWIRES;

};

#endif
