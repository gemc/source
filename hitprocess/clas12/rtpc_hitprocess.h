#ifndef RTPC_HITPROCESS_H
#define RTPC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

static const double PI=3.1415926535;

// constants to be used in the digitization routine
class rtpcConstants
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
	  // drift parameters 
	double z0[7], z2[7], z4[7];
	  // gain balance (example) keep for future or reference
	//double gain_balance[96][180]; // 180 rows, 96 cols
	
};


// Class definition
/// \class rtpc_HitProcess

class rtpc_HitProcess : public HitProcess
{
public:
	
	~rtpc_HitProcess(){;}
	
	// constants initialized with initWithRunNumber
	static rtpcConstants rtpcc;
	
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
	static HitProcess *createHitClass() {return new rtpc_HitProcess;}
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
private:

	// RTPC geometry parameters (constant) These should be consistant with PadVector.java in coatjava.
	const float PAD_W = 2.79;
	const float PAD_L = 4.0;
	const float PAD_S = 79.0; // old value was 80.0
	const float RTPC_L = 384.0;
	const float phi_per_pad = (2.0*PI)/180;
	
	// parameters for drift and diffustion equations for drift time, 
	// drift angle, and diffusion in z
	float a_t, b_t, c_t, diff_at, diff_bt;
	float a_phi, b_phi, c_phi, diff_aphi, diff_bphi;
	float a_z, b_z;
	
	// variables for storing drift times and diffusion in time
	//float sigma_t_2GEM2, sigma_t_2GEM3, sigma_t_2PAD, sigma_t_gap;
	const float sigma_t_2GEM2 = 8.72728;
	const float sigma_t_2GEM3 = 5.62223;
	const float sigma_t_2PAD = 7.58056;
	float sigma_t_gap;
	
	// variables for storing drift angle and diffusion in phi
	float phi_2GEM2 = 0.0; // 0.0416925;
	float phi_2GEM3 = 0.0; // 0.0416574;
	float phi_2PAD = 0.0; // 0.057566;
	float phi_2END;
	//float sigma_phi_2GEM2, sigma_phi_2GEM3, sigma_phi_2PAD, sigma_phi_gap;
	const float sigma_phi_2GEM2 = 0.00384579;
	const float sigma_phi_2GEM3 = 0.00160235;
	const float sigma_phi_2PAD = 0.00238653;
	float sigma_phi_gap;

	
	float z_cm;
	float TPC_TZERO; // What's this?
	
	map<int, double> timeShift_map;
	double shift_t;
	
};

#endif
