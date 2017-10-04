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
    string variation;
    string date;
    string connection;
    
    double driftVelocity[6];
    double miniStagger[6];
    double dcThreshold;
    int NWIRES;
    double dLayer[6];                              // ~cell size in each superlayer - one of Mac's core parameters
    
    // efficiency parameters for each superlayer
    double P1[6], P2[6], P3[6], P4[6], iScale[6];
    
    // smearing parameters for each sector / superlayer
    double smearP1[6][6], smearP2[6][6], smearP3[6][6], smearP4[6][6], smearScale[6][6];
    
    //  smearing parametes for the (random) time walk contributions:
    double smear_time_walk[3];
    
    //	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
    double vpar[4];
    
    // translation table
    TranslationTable TT;
    
    
    //parameters for time to distance:
    double deltanm[6][6], v0[6][6], delta_bfield_coefficient[6][6],tmaxsuperlayer[6][6];
    double deltatime_bfield_par1[6][6], deltatime_bfield_par2[6][6], deltatime_bfield_par3[6][6], deltatime_bfield_par4[6][6];
    double dmaxsuperlayer[6];
    
    
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
    
    // returns a time given a distance
    double calc_Time(double x, double dmax, double tmax, double alpha, double bfield, int sector, int superlayer);
    
    // returns time walks according to ionisation process:
    double time_walk_core(double x, double dmax, double epsilon, double R, double kappa, double v0);
    
    // returns random time walks:
    double time_rnd_core(double x, double f, double v0);
    
    G4ThreeVector psmear(G4ThreeVector p);
    
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
