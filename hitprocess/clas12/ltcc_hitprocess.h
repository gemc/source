#ifndef LTCC_HITPROCESS_H
#define LTCC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// constants to be used in the digitization routine
class ltccConstants
{
public:
	
	// database
	int    runNo;
	string date;
	string connection;
	char   database[80];
	
	// translation table
	TranslationTable TT;
	
	// add constants here
	double speMean[6][2][18] = {};
	double speSigma[6][2][18] = {};
	double timeOffset[6][2][18] = {};
	double timeRes[6][2][18] = {};

    // tdcs conversion factor
    double tdc_conv[6][2][18] = {};

};



// Class definition
class ltcc_HitProcess : public HitProcess
{
public:
	
	~ltcc_HitProcess(){;}
	
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
	static HitProcess *createHitClass() {return new ltcc_HitProcess;}
	
private:
	
	// constants initialized with initWithRunNumber
	static ltccConstants ltccc;
	
	void initWithRunNumber(int runno);
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
	double fadc_precision = 0.0625;  // 62 picoseconds resolution
	double convert_to_precision(double time) {
		return (int( time / fadc_precision ) * fadc_precision);
	}
	

};

#endif
