#ifndef FTOF_HITPROCESS_H
#define FTOF_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"


class ftofConstants
{
	public:

		// database
		int    runNo;	
		string variation;
		string date;
		string connection;
		char database[80];

	
		// There are 3 FTOF panels, 6 sectors so the constants are organized in
		// 6 + 3 dimensional arrays + 2 for Left and Right
	
		// number of paddles in each panel, in order p1a, p1b, p2
		int npaddles[3];
	
		// status:
		//	0 - fully functioning
		//	1 - noADC
		//	2 - noTDC
		//	3 - noADC, noTDC(PMTisdead)
		//  5 - any other reconstruction problem
		vector<int> status[6][3][2];
	
		// effective velocity
		vector<double> veff[6][3][2];
	
		// attenuation length comes from a linear parameterization
		// based on the counter length in cm
		// depends on the panel
		vector<double> attlen[6][3][2];
	
		// de/dx = 2MeV / g/ cm3 for MIP in the FTOF scintillators
		double dEdxMIP;
	
		// minimum ionizing calibration peak
		vector<double> countsForMIP[6][3][2];
	
		// time walk correction is parameterized with two coefficients
		vector<double> twlk[6][3][6];
	
		// time resolution parameterized as sigma0^2 + sigma1^2/N
		// where N is number of photoelectrons reaching the PTMS
		double sigma0[3], sigma1[3];
		double nphePerMevReachingPMT;
	
};


// Class definition
/// \class ftof_HitProcess
/// <b> Forward Time of Flight Hit Process Routine</b>\n\n
/// The Calibration Constants are:\n
/// - VEF is the effective velocity of propogation in the scintillator
class ftof_HitProcess : public HitProcess
{
	public:
	
		~ftof_HitProcess(){;}
	
		// constants initialized with initWithRunNumber
		static ftofConstants ftc;
	
		void initWithRunNumber(int runno);

		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);

		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
		// creates the HitProcess
		static HitProcess *createHitClass() {return new ftof_HitProcess;}
	
	
};



#endif




