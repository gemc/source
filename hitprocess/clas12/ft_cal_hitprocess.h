#ifndef ft_cal_HITPROCESS_H
#define ft_cal_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"


// Class definition
class ft_cal_HitProcess : public HitProcess
{
public:

	~ft_cal_HitProcess(){;}

	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);

	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);

	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);

	// creates the HitProcess
	static HitProcess *createHitClass() {return new ft_cal_HitProcess;}

	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

};

#endif
