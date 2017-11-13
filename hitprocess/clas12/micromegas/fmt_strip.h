#ifndef fmt_strip_H
#define fmt_strip_H 1

// geant4
#include "G4ThreeVector.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;
#include <string>
#include <vector>
using namespace std;

#include "Lorentz.h"

class fmtConstants
{
public:

	// database
	int    runNo;
	string variation;
	string date;
	string connection;
	char   database[80];
	double SigmaDrift;	 // Max transverse diffusion value (GEMC value)
	double hDrift;			 // Size of the drift gap
	double hStrip2Det;	 // Distance between strips and the middle of the conversion gap (~half the drift gap)
	int nb_sigma;            // Define the number of strips to look around the closest strip
	double ThetaL;		    // Lorentz angle
	double Theta_Ls;         //Angle of the Lorentz drift in xy plane 
	double fieldScale; 	 // Scaling of the field - 1 means a 5 Tesla field along z

	// THE GEOMETRY CONSTANTS
	const static int NLAYERS = 6  ;	// 6 DISKS

	// Detector geometrical characteristics
	vector<double> Z0;       // z of the upstream part of the layer
	vector<double> alpha;    // strip angles of layers
	double R_max;            // outer radius of strip part
	double R_min;            // inner radius of strip part
	double R_IR=84.05;             //Radius of the inner region of the strips
	double pitch;   // the width of the corresponding group of strips
	int N_str;               // Number of strips for 1 layer
	int N_halfstr;           // number of bottom strips in the central part
	int N_sidestr;           // number of strips one side
	double y_central;        // Limit the central part for stip finding

	// Detector response characteristics
	double HV_DRIFT[NLAYERS]; //Need to know the HV to compute the lorentz angle...0 upstream
	double HV_STRIPS_IN[NLAYERS]; //Might be needed for the gain
	double HV_STRIPS_OUT[NLAYERS]; //Might be needed for the gain
	vector<vector<int> >     STRIP_STATUS; // Say if strip is dead or alive
	vector<vector<double> >     STRIP_GAIN; // Give the gain where the strip is
	vector<vector<double> >     STRIP_EFFICIENCY; // Give the efficiency (correlated to gain fluctuation)

	double w_i=20; //ionization potential assumed to be 25 eV
	Lorentz Lor_Angle;

};

class fmt_strip
{
public:
	double sigma_td;         // current value of the transverse diffusion
	int nb_sigma;            // Number of strips to study around the closest one
	int Nel;                 // number of electrons (Nt) for a given hit
	double y_real; 		 // y position of the hit after transverse diffusion
	double x_real;		 // x position of the hit after transverse diffusion
	double strip_x;          // strip_x is the middle of the strips
	double strip_y;          // strip_y is the position of the strips
	double strip_length;     // length of the strip

	vector<double> FindStrip( int layer, int sector, double x, double y, double z, double Edep, fmtConstants fmtc);   // Strip Finding Routine
	void Carac_strip(int strip, fmtConstants fmtc); //length of the strip
	double Weight_td(int strip, double x, double y, double z, fmtConstants fmtc); //Compute the fraction of Nel falling onto the strip, depending on x,y in the FMT coordinate system
	double GetBinomial(double n, double p);//CLHEP Binomial has a weird limit condition which returns -1 instead of 0 when n*p=0
 };

#endif
