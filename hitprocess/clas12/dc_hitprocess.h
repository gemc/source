#ifndef DC_HITPROCESS_H
#define DC_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// constants to be used in the digitization routine
// warning: since NWIRES and ministagger are also used by processID, the plugin loading
// has to happen before the first event is processed. In other words,
// initializeDCConstants(1) - or remove the 	if(runno == -1) return dcc;
class dcConstants
{
public:
	
	// database
	int    runNo;
	string connection;
	char   database[80];
	
	int    fieldPolarity;
	//	double driftVelocity[6];
	double miniStagger[6];
	double dcThreshold;
	int NWIRES = 112;
	int NLAYERS= 6;
	double dLayer[6];                              // ~cell size in each superlayer - one of Mac's core parameters
	
	// efficiency parameters for each superlayer
	double P1[6][6], P2[6][6], P3[6][6], P4[6][6], iScale[6][6];
	
	// smearing parameters for each sector / superlayer
	double smearP0[6][6], smearP1[6][6], smearP2[6][6], smearP3[6][6], smearP4[6][6];
	
	//	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
	double vpar[4];
	
	// translation table
	TranslationTable TT;
	
	
	//parameters for time to distance:
	double deltanm[6][6], v0[6][6], delta_bfield_coefficient[6][6],tmaxsuperlayer[6][6];
	double deltatime_bfield_par1[6][6], deltatime_bfield_par2[6][6], deltatime_bfield_par3[6][6], deltatime_bfield_par4[6][6];
	double vmid[6][6], R[6][6];
	double dmaxsuperlayer[6];
	
	// sector, SL, slot, cable
	double T0Correction[6][6][7][6];
	
	
	double get_T0(int sectorI, int superlayerI, int layerI, int nwire) {
		
		int slot = ((nwire - 1) / 16) + 1;
		
		int wire1to16I = ((nwire - 1) % 16);
		int cable = CableID[layerI][wire1to16I];
		
		double t0corr = T0Correction[sectorI][superlayerI][slot - 1][cable - 1];
		
		//		cout << "  sectorI: " << sectorI << ", superlayerI: " << superlayerI << ",  layerI: " << layerI << ",  nwire: "
		//		<< nwire << ", slot: " << slot << ", cable: " << cable << "  T0: " << t0corr << endl;
		
		return t0corr;
	}
	
	
	int CableID[6][16] = {
		//[nLayer][nLocWire] => nLocWire=16, 7 groups of 16 wires in each layer
		{1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6}, //Layer 1
		{1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6}, //Layer 2
		{1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6}, //Layer 3
		{1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6}, //Layer 4
		{1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6}, //Layer 5
		{1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6}, //Layer 6
	};
	//===> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
	// (Local wire ID: 0 for 1st, 16th, 32th, 48th, 64th, 80th, 96th wires)

	
	
};






// Class definition
class dc_HitProcess : public HitProcess
{
public:
	
	~dc_HitProcess(){;}
	
	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);
	
	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);
	
	
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
	// creates the HitProcess
	static HitProcess *createHitClass() {return new dc_HitProcess;}
	
	// returns a time given a distance: old exponential function
	double calc_Time_exp(double x, double dmax, double tmax, double alpha, double bfield, int sector, int superlayer);
	
	// returns a time given a distance: neew polynomial function
	double calc_Time(double x, double dmax, double tmax, double alpha, double bfield, int sector, int superlayer);
	
	// returns time walks according to ionisation process:
	double doca_smearing(double x, double beta, int sector, int superlayer);
	
	G4ThreeVector psmear(G4ThreeVector p);

        G4ThreeVector wireLxyz(int layer, int wire, double dlayer, double dwire);

        double doca(G4ThreeVector pos, int layer, int wire, double dlayer, double dwire);

private:
	
	// constants initialized with initWithRunNumber
	static dcConstants dcc;
	
	void initWithRunNumber(int runno);
	
	//Is not used in --> Will be removed in gemc 3.0
	//********************************************************************************************
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
	
	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);
	
	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	//********************************************************************************************
	
};

#endif
