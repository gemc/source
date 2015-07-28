#ifndef EC_HITPROCESS_H
#define EC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// ec constants
// these are loaded with initWithRunNumber
class ecConstants
{
	public:
		// runNo is mandatory variable to keep track of run number changes
		int runNo;
	
		double NSTRIPS;              // Number of strips
		double attlen;               // Attenuation Length (mm)
		double TDC_time_to_channel;  // conversion from time (ns) to TDC channels.
		double ECfactor;             // number of p.e. divided by the energy deposited in MeV; value taken from gsim. see EC NIM paper table 1.
		int TDC_MAX;                 // max value for EC tdc.
		double ec_MeV_to_channel;    // conversion from energy (MeV) to ADC channels
};


// Class definition
class ec_HitProcess : public HitProcess
{
	public:

		~ec_HitProcess(){;}

		// constants initialized with initWithRunNumber
		static ecConstants ecc;
	
		void initWithRunNumber(int runno);

		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);
		
		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
		// creates the HitProcess
		static HitProcess *createHitClass() {return new ec_HitProcess;}
	
};





#endif
