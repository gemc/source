#ifndef CND_HITPROCESS_H
#define CND_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition
class cnd_HitProcess : public HitProcess
{
	public:
	
		~cnd_HitProcess(){;}
	
		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);
		
		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);

		// creates the HitProcess
		static HitProcess *createHitClass() {return new cnd_HitProcess;}

		double BirksAttenuation(double,double,int,double);
};

#endif
