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
		string variation;
	        string date;
		string connection;
		char   database[80];
		
		double NSTRIPS;                 // Number of strips
		double attlen[3][36][2][3][6];  // Attenuation Length (mm)
		double attl;
		double TDC_time_to_evio;        // Conversion from time (ns) to EVIO TDC format
		double ADC_MeV_to_evio;         // Conversion from energy (MeV) to EVIO FADC250 format
		double PE_yld;                  // Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
		double veff;                    // Effective velocity of scintillator light (mm/ns)
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
