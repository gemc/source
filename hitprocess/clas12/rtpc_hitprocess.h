#ifndef RTPC_HITPROCESS_H
#define RTPC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition
/// \class rtpc_HitProcess

class rtpc_HitProcess : public HitProcess
{
public:
	
	~rtpc_HitProcess(){;}
	
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
	
public:
	// RTPC geometry parameters
	float PAD_W, PAD_L, PAD_S, RTPC_L;
	float phi_per_pad;
	
	// parameters for drift and diffustion equations for drift time, 
	// drift angle, and diffusion in z
	float a_t, b_t, c_t, diff_at, diff_bt;
	float a_phi, b_phi, c_phi, diff_aphi, diff_bphi;
	float a_z, b_z;
	
	// variables for storing drift times and diffusion in time
	float t_2GEM2, t_2GEM3, t_2PAD, t_2END;
	float sigma_t_2GEM2, sigma_t_2GEM3, sigma_t_2PAD, sigma_t_gap;
	
	// variables for storing drift angle and diffusion in phi
	float phi_2GEM2, phi_2GEM3, phi_2PAD, phi_2END;
	float sigma_phi_2GEM2, sigma_phi_2GEM3, sigma_phi_2PAD, sigma_phi_gap;
	
	float z_cm;
	float TPC_TZERO;
	
	map<int, double> timeShift_map;
	double shift_t;
	
};

#endif
