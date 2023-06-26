#ifndef RICH_HITPROCESS_H
#define RICH_HITPROCESS_H 1

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>

// gemc headers
#include "HitProcess.h"

// geant4 headers
#include "Randomize.hh"
#include "G4Poisson.hh"

/* From Marco: class to generate the single photon signal of MAPMT and readout electronics 
   8500 -> Hamamatsu H8500
   12700 -> Hamamatsu H12700
   12701 -> Hamamatsu H12700 but with old (and wrong) stage gain calculation 
*/


class RichPixel {

private:

  /* charge of the electron in fC */
  const double Qe = 1.602e-4;
  
  //TRandom3 rnd;
  HepRandom rnd;
  /* ----------------------- */
  /* MAPMT parameters */
  int PmtType = 0;
  double Gain = 3e6;
  double G1 = 1;
  double GN = 1;
  double d1Ratio = 1;
  int nStages = 1;

  /* ----------------------- */
  /* MAROC parameters */

  /* MAROC charge saturation (fC) */
  const double MarocMaxQ = 2500;

  
  /* ADC conversion factors */

  /* conversion from fC to MAROC binary ADC units */
  const double DAC = 0.82;

  /* Pedestal in DAC units */
  const int Pedestal = 0;



  /* TDC conversion factors */

  /* scaling factor to convert the threshold used for the timing to the one used to cut hits */
  double ThresholdScale = 1;

  /* from charge to time */
  double q0_t = 0;
  double p_t[5] = {0, 0, 0, 0, 0};
  double m_t = 0;
  double q_t =0;

  /* from charge to duration */
  double p0_d = 0;
  double p1_d = 0;

  double alphaD = 0.;

  /* MAROC settings */

  /* Default threshold at +25 DAC */
  double StdTrhesholdDAC = 25;

  /* MAROC shaping gain */
  double MarocG = 1;

  /* Relative Threshold setting (between 1 and 2) */
  double MarocThrF = 1;

  /* TDC threshold in charge */
  double MarocThrCharge = 0;

  /* Time offset */
  double TimeOffset = 0;

  /* Time sigma resolution for the gaussian smearing */
  double TimeResol = 1;

  /* Output Signal quantities */
  int npe;

  double qadc;
  int ADC;
  
  double qtdc;
  double start_time;
  double true_t1;
  double t1, t2;
  double duration;


  void GenerateNpe(int n0);

  void ChargeToTime();
  void ChargeToDuration();


public:

  RichPixel(int type=12700);

  void InitReadout(int pmttype, double pmtgain, double marocgain, double marocthr);
  void InitPmt(int type=12700, double g=3e6);
  void InitMaroc(double g=1, double th=1);


  void Clear();

  int get_PmtType() { return PmtType; }
  int get_Npe() { return npe; }
  double get_ChargeADC() { return qadc; }
  int get_ADC() { return ADC; }

  double get_ChargeTDC() { return qtdc; }
  double get_TrueT1() { return true_t1; }
  int get_T1() { return (int)t1; }
  int get_T2() { return (int)t2; }
  int get_Duration() { return (int)(t2-t1); }


  int GenerateADC(int n0);
  bool GenerateTDC(int n0, double t0);
  bool GenerateTDCfromADC(double qadc, double t0);

  void PrintPmt();
  void PrintMaroc();
};

// constants to be used in the digitization routine
class richConstants
{
public:
	
        // database                                                                                                                                              
        int    runNo;
        string date;
        string connection;
        char   database[80];
        string variation;
        RichPixel *richPixel;
        // translation table                                                                                                                                     
        TranslationTable TT;

  // add constants here                                                                                                                                          
  const static int npmt = 391;
  const static int npixel = 64;
  double timewalkCorr_D0[npmt];
  double timewalkCorr_m1[npmt];
  double timewalkCorr_m2[npmt];
  double timewalkCorr_T0[npmt];
  double timeOffsetCorr[npmt];

	
};



// Class definition
class rich_HitProcess : public HitProcess
{
public:
	
	~rich_HitProcess(){;}
	
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
	static HitProcess *createHitClass() {return new rich_HitProcess;}
	
private:
	
	// constants initialized with initWithRunNumber
	static richConstants richc;
	
	void initWithRunNumber(int runno);
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
};

#endif
