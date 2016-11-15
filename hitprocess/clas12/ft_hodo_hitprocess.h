#ifndef ft_hodo_HITPROCESS_H
#define ft_hodo_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// constants to be used in the digitization routine
class ftHodoConstants
{
public:

	// database
	int    runNo;
	string variation;
	string date;
	string connection;

	// translation table
	TranslationTable TT;

	// add constants here

};



// Class definition
class ft_hodo_HitProcess : public HitProcess
{
public:

	~ft_hodo_HitProcess(){;}

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
	static HitProcess *createHitClass() {return new ft_hodo_HitProcess;}

private:

	// constants initialized with initWithRunNumber
	static ftHodoConstants fthc;

	void initWithRunNumber(int runno);

	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

};

#endif
