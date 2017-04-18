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
    
    double sigma_0      ;  // transverse diffusion
    double w_i          ;  // ionization energy
    
    double rmin ;          // inner radius of disks
    double rmax ;          // outer radius of disks
    double pitch;          // strip pitch
    int    nstrips;        // number of strips per layer
	int    nb_sigma;       // To define the number of strips to look at around the closest one
    //	voltage signal parameters, using double gaussian + delay (function DGauss, need documentation for it)
    double vpar[4];
};


class ftm_strip
{
public:
	vector<double> FindStrip( int layer, double x, double y, double z, double Edep, detector Detector,ftmConstants ftmcc);   // Strip Finding Routine
    int    get_strip_ID(double x, double y,ftmConstants ftmcc);   // return strip number based on coordinate
    double get_strip_X(double x, double y,ftmConstants ftmcc);    // return strip x coordinate
//    double get_strip_L(double x, double y,ftmConstants ftmcc);    // return strip length
};

#endif