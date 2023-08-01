#ifndef band_HITPROCESS_H
#define band_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// constants to be used in the digitization routine
class bandHitConstants
{
public:
	
	// database
	int    runNo;
	string date;
	string connection;
	string variation;
	char   database[80];
	
	int nsector;
	int nlayer;
	int ncomp;
	
	double mev_adc[6][6][7];
	double eff_vel_tdc[6][6][7];
	double eff_vel_fadc[6][6][7];
	double atten_len[6][6][7];
	double tdc_offset[6][6][7];
	double tdc_resolution[6][6][7];
	
	//	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
	double vpar[4];

    // tdc conversion factor
    double tdcconv;


};



// Class definition
class band_HitProcess : public HitProcess
{
public:
	
	~band_HitProcess(){;}
	
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
	static HitProcess *createHitClass() {return new band_HitProcess;}
	
	double MeVtoMeVee(int PID, int Z, double E_MeV );
	double BirksAttenuation(double destep, double stepl, int charge, double birks);
	
private:
	
	// constants initialized with initWithRunNumber
	static bandHitConstants bhc;
	
	void initWithRunNumber(int runno);
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

    double fadc_precision = 0.0625;  // 62 picoseconds resolution
    int convert_to_precision(double time) {
        return (int( time / fadc_precision ) / fadc_precision);
    }
};

#endif
