#ifndef MYATOF_HITPROCESS_H
#define MYATOF_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"


class myatofConstants
{
public:
	
	// Database parameters
	int    runNo;
	string date;
	string connection;
	char   database[80];

	// translation table
	TranslationTable TT;
};


// Class definition
/// \class atof_HitProcess
/// <b> Alert Time of Flight Hit Process Routine</b>\n\n

class myatof_HitProcess : public HitProcess
{
public:
	
	~myatof_HitProcess(){;}
	
	// constants initialized with initWithRunNumber
	static myatofConstants atc;
	
	void initWithRunNumber(int runno);
	
	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);
	
	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);
	
	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);
	
	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
	// creates the HitProcess
	static HitProcess *createHitClass() {return new myatof_HitProcess;}
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
};



#endif




