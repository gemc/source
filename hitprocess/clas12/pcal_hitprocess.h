#ifndef PCAL_HITPROCESS_H
#define PCAL_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// pc constants
// these are loaded with initWithRunNumber
class pcConstants
{
public:

	// runNo is mandatory variable to keep track of run number changes
	int    runNo;
	string variation;
	string date;
	string connection;
	char   database[80];

        static const int nsect  = 6;  // Number of sectors
        static const int nlayer = 9;  // layer=1-3 (PCAL) 4-6 (ECinner) 7-9 (ECouter)
        static const int nview  = 3;  // Number of views, U,V and W
	
	// For strip dependent constants read from CCDB
	// Array [6][9][3] -> sector,layer,view sector=1-6 layer=1-3 (PCAL) 4-6 (ECinner) 7-9 (ECouter) view=1-3 (U,V,W)

	//attlen: attenuation length
	vector<double> attlen[nsect][nlayer][nview];
	
        //gain: pmt gain
        vector<double> gain[nsect][nlayer];
	
	//timing: TDC calibration constants
	vector<double> timing[nsect][nlayer][5];

	
        // ======== FADC Pedestals and sigmas ===========
        double pedestal[nsect][nlayer][nview] = {};
	double pedestal_sigm[nsect][nlayer][nview] = {};
	
	double TDC_time_to_evio;     // Conversion from time (ns) to EVIO TDC format
	double ADC_GeV_to_evio;      // Conversion from energy (MeV) to EVIO FADC250 format
	double veff;                 // Effective velocity of scintillator light (mm/ns)
	double pmtPEYld;             // Number of p.e. divided by the energy deposited in MeV.
	double pmtQE;
	double pmtDynodeGain;
	double pmtDynodeK;
	double pmtFactor;

	// voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
	double vpar[4];

	// translation table
	TranslationTable TT;
};

// Class definition
class pcal_HitProcess : public HitProcess
{
public:

	~pcal_HitProcess(){;}

	// constants initialized with initWithRunNumber
	static pcConstants pcc;

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
	static HitProcess *createHitClass() {return new pcal_HitProcess;}

	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
};

#endif
