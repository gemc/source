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
    char   database[80];

    // translation table
    TranslationTable TT;

	// add constants here
    // For crystal dependent constants read from CCDB
    // Array [8][2] -> sector,layer,component
    
    // status:
    //	0 - fully functioning
    //	1 - noisy channel
    //	3 - dead channel
    //  5 - any other issue
    vector<int> status[8][2];
    
    // noise
    vector<double>  pedestal[8][2];
    vector<double>  pedestal_rms[8][2];
    vector<double>  gain_pc[8][2];
    vector<double>  gain_mv[8][2];
    vector<double>  npe_threshold[8][2];
    
    // energy
    vector<double>  mips_charge[8][2];
    vector<double>  mips_energy[8][2];
    vector<double>  preamp_gain[8][2];
    
    // time
    vector<double>  time_offset[8][2];
    vector<double>  time_rms[8][2];
    
    // fadc parameters
    double  ns_per_sample;
    double  fadc_input_impedence;
    double  fadc_LSB;
    double  time_to_tdc;
    double  tdc_max;
        
    //	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
    double vpar[4];

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
