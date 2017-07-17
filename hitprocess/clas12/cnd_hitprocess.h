#ifndef CND_HITPROCESS_H
#define CND_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"


// constants to be used in the digitization routine
class cndConstants
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
	/* vector<int> status[24][3][2]; */
	/* vector<double> veff[24][3][2]; */
	/* vector<double> att_length[24][3][2]; */
	/* vector<double> time_offset_LR[24][3][2]; */
	/* vector<double> time_offset_layer[24][3][2]; */
	/* vector<double> uturn_t[24][3][2]; */
	/* vector<double> uturn_e[24][3][2]; */
	/* vector<double> ecal[24][3][4]; */

	/*	int status[24][3][2];
	double slope[24][3][2];
	double veff[24][3][2];
	double att_length[24][3][2];
	double time_offset_LR[24][3][1];
	double time_offset_layer[24][3][1];
	double uturn_t[24][3][1];
	double uturn_e[24][3][1];
	double ecalD[24][3][2];
	double ecalN[24][3][2];
	*/

	int status_L[24][3][2];                                                                       int status_R[24][3][2];             
        double slope_L[24][3][2];                                                                     double slope_R[24][3][2];               
        double veff_L[24][3][2];                                                                      double veff_R[24][3][2]; 
	double attlen_L[24][3][2]; 
        double attlen_R[24][3][2];                                                            
        double time_offset_LR[24][3][1];                                                              double time_offset_layer[24][3][1];                                                           double uturn_tloss[24][3][1];                                                                 double uturn_e[24][3][1];                                                                     double mip_dir_L[24][3][2];                                                                   double mip_dir_R[24][3][2];          
        double mip_indir_L[24][3][2];                                                         
        double mip_indir_R[24][3][2];

};



// Class definition
class cnd_HitProcess : public HitProcess
{
public:

	~cnd_HitProcess(){;}

	static cndConstants cndc;

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
	static HitProcess *createHitClass() {return new cnd_HitProcess;}

	double BirksAttenuation(double,double,int,double);

private:


	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
};

#endif
