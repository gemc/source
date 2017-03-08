#ifndef ftm_strip_H
#define ftm_strip_H 1

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

#include <vector>
#include <string>
using namespace std;

// gemc headers
#include "HitProcess.h"

// ftm constants
// these are loaded with initWithRunNumber
class ftmConstants
{
public:
    
    // runNo is mandatory variable to keep track of run number changes
    int runNo;
    string variation;
    string date;
    string connection;
    char   database[80];
    
    double sigma_td_max ;  // // very small transverse diffusion (temporary)
    double w_i          ;  // ionization energy
    
    double rmin ;          // inner radius of disks
    double rmax ;          // outer radius of disks
    double pitch;          // strip pitch
    int    nstrips;        // number of strips per layer
    
    //	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
    double vpar[4];
};


class ftm_strip
{
public:
/*
	double pitch;
	double Pi;

	double interlayer;       // distance between 2 layers of a superlayer
	double intersuperlayer;  // distance between 2 superlayers
	double Rmin;             // inner radius of disks
	double Rmax;             // outer radius of disks
//	double Z_1stlayer;       // z position of the 1st layer

	vector<double> Z0;       // z of the upstream part of the layer
	vector<double> R;        // radii of layers
	vector<double> MidTile;  // mid angle of the sector
	int Nstrips;             // Number of strips for 1 card
//	double hDrift;           // Size of the drift gap
//	double hStrip2Det;       // distance between strips and the middle of the conversion gap (~half the drift gap)
	double sigma_td_max;     // maximum value of the transverse diffusion
	double sigma_td;         // current value of the transverse diffusion
	double w_i;              // mean ionization potential
	int Nel;                 // number of electrons (Nt) for a given hit
	double x_real, y_real; // x,y of the hit after transverse diffusion
 */
//	void fill_infos(detector);

	vector<double> FindStrip( int layer, double x, double y, double z, double Edep, detector Detector,ftmConstants ftmcc);   // Strip Finding Routine

};

#endif

