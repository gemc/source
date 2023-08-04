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

public:

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
  // anode->maroc and pmt->board constants

  // anode->maroc
  int anodeToMaroc[64] = {32, 30, 31, 29, 33, 35, 34, 36, 28, 26, 27, 25, 37,
       39, 38, 40, 24, 22, 23, 21, 41, 43, 42, 44, 20, 18,
       19, 17, 45, 47, 46, 48, 16, 14, 15, 13, 49, 51, 50,
       52, 12, 10, 11,  9, 53, 55, 54, 56,  8,  6,  7,  5,
       57, 59, 58, 60,  4,  2,  3,  1, 61, 63, 62, 64};
  // pmt->tile
  int pmtToTile[397] = {1,   1,   1,   2,   2,   2,   3,   3,   4,   4,   4,
		     5,   5,   6,   6,   6,   7,   7,   8,   8,   8,   9,
		     9,   9,  10,  10,  10,  11,  11,  11,  12,  12,  13,
		     13,  13,  14,  14,  14,  15,  15,  16,  16,  16,  17,
		     17,  17,  18,  18,  18,  19,  19,  20,  20,  20,  21,
		     21,  21,  22,  22,  22,  23,  23,  23,  24,  24,  25,
		     25,  25,  26,  26,  26,  27,  27,  27,  28,  28,  29,
		     29,  29,  30,  30,  30,  31,  31,  31,  32,  32,  32,
		     33,  33,  34,  34,  34,  35,  35,  35,  36,  36,  36,
		     37,  37,  37,  38,  38,  38,  39,  39,  40,  40,  40,
		     41,  41,  41,  42,  42,  42,  43,  43,  43,  44,  44,
		     45,  45,  45,  46,  46,  46,  47,  47,  47,  48,  48,
		     48,  49,  49,  49,  50,  50,  51,  51,  51,  52,  52,
		     52,  53,  53,  53,  54,  54,  54,  55,  55,  55,  56,
		     56,  56,  57,  57,  58,  58,  58,  59,  59,  59,  60,
		     60,  60,  61,  61,  61,  62,  62,  62,  63,  63,  64,
		     64,  64,  65,  65,  65,  66,  66,  66,  67,  67,  67,
		     68,  68,  68,  69,  69,  69,  70,  70,  71,  71,  71,
		     72,  72,  72,  73,  73,  73,  74,  74,  74,  75,  75,
		     75,  76,  76,  76,  77,  77,  77,  78,  78,  79,  79,
		     79,  80,  80,  80,  81,  81,  81,  82,  82,  82,  83,
		     83,  83,  84,  84,  84,  85,  85,  86,  86,  86,  87,
		     87,  87,  88,  88,  88,  89,  89,  89,  90,  90,  90,
		     91,  91,  91,  92,  92,  92,  93,  93,  94,  94,  94,
		     95,  95,  95,  96,  96,  96,  97,  97,  97,  98,  98,
		     98,  99,  99,  99, 100, 100, 100, 101, 101, 101, 102,
		     102, 103, 103, 103, 104, 104, 104, 105, 105, 105, 106,
		     106, 106, 107, 107, 107, 108, 108, 108, 109, 109, 109,
		     110, 110, 111, 111, 111, 112, 112, 112, 113, 113, 113,
		     114, 114, 114, 115, 115, 115, 116, 116, 116, 117, 117,
		     117, 118, 118, 118, 119, 119, 120, 120, 120, 121, 121,
		     121, 122, 122, 122, 123, 123, 123, 124, 124, 124, 125,
		     125, 125, 126, 126, 126, 127, 127, 127, 128, 128, 128,
		     129, 129, 130, 130, 130, 131, 131, 131, 132, 132, 132,
		     133, 133, 133, 134, 134, 134, 135, 135, 135, 136, 136,
		     136, 137, 137, 137, 138, 138, 139, 139, 139, 140, 140,
		     140};
  // pmt->position on readout board (for maroc channel -> final component number)
  int pmtToTilePosition[397] = {1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1,
			     2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2,
			     3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
			     1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3,
			     1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3,
			     1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
			     1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
			     3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
			     3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2,
			     3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1,
			     2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1,
			     2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1,
			     2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 2, 3, 1, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
			     1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
			     1, 2, 3, 1, 2, 3};

  // quantum efficiency table
  const static int nQEbins=87;
  double Ene[nQEbins]={6.148,     5.387,
				      4.624,      4.558,      4.510,      4.432,      4.372,      4.285,      4.201,      4.081,      3.980,      3.815,
				      3.726,      3.591,      3.503,      3.410,      3.314,      3.214,      3.128,      3.054,      2.983,      2.902,
				      2.845,      2.783,      2.724,      2.668,      2.624,      2.572,      2.531,      2.487,      2.454,      2.431,
				      2.408,      2.381,      2.364,      2.342,      2.325,      2.305,      2.284,      2.260,      2.240,      2.221,
				      2.202,      2.180,      2.165,      2.147,      2.136,      2.119,      2.108,      2.091,      2.084,      2.071,
				      2.058,      2.048,      2.041,      2.032,      2.022,      2.013,      2.003,      1.997,      1.991,      1.982,
				      1.975,      1.969,      1.963,      1.957,      1.951,      1.946,      1.940,      1.934,      1.931,      1.922,
				      1.919,      1.914,      1.908,      1.902,      1.897,      1.894,      1.888,      1.883,      1.877,      1.872,
				      1.867,      1.861,      1.858,      1.853,      1.848};
  
  double QE[nQEbins]={0.000,     0.000,
			0.11241,    0.12738,    0.14272,    0.16358,    0.18120,    0.20301,    0.22234,    0.24351,    0.25775,    0.27282,
			0.27594,    0.27594,    0.27594,    0.27282,    0.26974,    0.26669,    0.26070,    0.25484,    0.24629,    0.23534,
			0.22745,    0.21489,    0.20301,    0.18963,    0.17713,    0.16358,    0.15279,    0.13951,    0.12738,    0.11500,
			0.10381,    0.09161,    0.08177,    0.07298,    0.06514,    0.05881,    0.05249,    0.04685,    0.04181,    0.03775,
			0.03369,    0.03007,    0.02684,    0.02396,    0.02163,    0.01887,    0.01703,    0.01503,    0.01342,    0.01184,
			0.01045,    0.00933,    0.00832,    0.00726,    0.00656,    0.00572,    0.00511,    0.00456,    0.00402,    0.00351,
			0.00317,    0.00276,    0.00247,    0.00215,    0.00192,    0.00169,    0.00150,    0.00130,    0.00118,    0.00103,
			0.00093,    0.00082,    0.00072,    0.00063,    0.00056,    0.00050,    0.00043,    0.00039,    0.00035,    0.00030,
			0.00027,    0.00023,    0.00021,    0.00018,    0.00016};

  // two types of pmts used in sector 4 rich
  int pmtType[391] = {12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 12700, 12700, 12700, 8500, 8500, 8500, 8500, 8500, 8500, 8500, 
		      8500, 8500, 8500, 8500, 8500, 8500, 8500, 8500, 8500, 12700, 12700, 12700, 8500, 8500, 12700, 12700, 12700, 12700, 
		      12700, 12700, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 12700, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 12700, 12700, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500,
		      8500, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 8500,
		      8500, 8500, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500,
		      12700, 12700, 12700, 8500, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 8500, 8500, 8500, 12700, 12700, 12700, 8500, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 12700, 
		      12700, 12700, 8500, 8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 8500, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500,
		      8500, 8500, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 8500, 8500, 12700,
		      12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700, 12700,
		      12700, 12700, 12700, 8500, 8500, 8500, 8500, 8500};
  
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

        // RICH specific functions 
        int getPixelNumber(G4ThreeVector  Lxyz);
        G4ThreeVector getPixelCenter(int pixel);        
};

#endif
