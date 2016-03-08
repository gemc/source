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
		int    runNo;
		string variation;
		string date;
		string connection;
		char   database[80];

		// For strip dependent constants read from CCDB
		// Array [6][9][3] -> sector,layer,view sector=1-6 layer=1-3 (PCAL) 4-6 (ECinner) 7-9 (ECouter) view=1-3 (U,V,W)

		//attlen: attenuation length
		vector<double> attlen[6][9][3];
		
		double NSTRIPS;             // Number of strips
		double TDC_time_to_evio;    // Conversion from time (ns) to EVIO TDC format
		double ADC_MeV_to_evio;     // Conversion from energy (MeV) to EVIO FADC250 format
		double veff;                // Effective velocity of scintillator light (mm/ns)
		double pmtPEYld;            // Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
		double pmtQE;               // Quantum efficiency of PMT
		double pmtDynodeGain;       // PMT dynode gain
		double pmtDynodeK;          // PMT dynode secondary emission statistics factor: K=0 (Poisson) K=1 (exponential) 
		double pmtFactor;           // Contribution to FWHM from PMT statistical fluctuations.		  
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
