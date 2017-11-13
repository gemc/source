#ifndef bmt_strip_H
#define bmt_strip_H 1

// geant4
#include "G4ThreeVector.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// c++
#include <vector>
#include <string>
#include "Lorentz.h"
using namespace std;

/// \class bmt_strip
/// <b> bmt_strip </b>\n\n
/// Micromegas strip finding routine\n

// constants to be used in the digitization routine
// warning: since fieldScale is also used by processID, the plugin loading
// has to happen before the first event is processed.
class bmtConstants
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
	double Theta_Ls_Z;         //Angle between the Lorentz deviation and the strip.
	double Theta_Ls_C;         //Angle between the Lorentz deviation and the strip.
	//double fieldScale; 	 // Scaling of the field - 1 means a 5 Tesla field along z

	// THE GEOMETRY CONSTANTS
	const static int NLAYERS = 6  ;	// 6 layers of MM
	const static int NSECTORS = 3  ;	// 3 tiles per layer of MM

	// Detector geometrical characteristics
	double RADIUS[NLAYERS] ; 	// the radius of the Z detector in mm
	int AXIS[NLAYERS] ;             // 0 if C-detector, 1 if Z detector
	int NSTRIPS[NLAYERS] ; 	       // the number of strips
	double LENGTH[NLAYERS] ; 	       // the strip length in mm
	double ZMIN[NLAYERS] ; 	       // PCB upstream extremity mm
	double ZMAX[NLAYERS] ; 	       // PCB downstream extremity mm
	double EDGE1[NLAYERS][NSECTORS]; // the angle of the first edge of each PCB detector A, B, C
	double EDGE2[NLAYERS][NSECTORS]; // the angle of the second edge of each PCB detector A, B, C
	
	vector<vector<int> >     GROUP;   // Number of strips with same width
	vector<vector<double> >  PITCH;   // the width of the corresponding group of strips

	// Detector response characteristics
	double HV_DRIFT[NLAYERS][NSECTORS]; //Need to know the HV to compute the lorentz angle
	double HV_STRIPS[NLAYERS][NSECTORS]; //Might be needed for the gain
	vector<vector<int> >     STRIP_STATUS; // Say if strip is dead or alive
	vector<vector<double> >     STRIP_GAIN; // Give the gain where the strip is
	vector<vector<double> >     STRIP_EFFICIENCY; // Give the efficiency (correlated to gain fluctuation)

	double w_i=20; //ionization potential assumed to be 25 eV... 
	double density=1.86e-3; //g*cm-3
	double np=30.27; //cm-1 number of primary ionization
	double nt=104.1; //cm-1 number of total ionizations
	double averA=36.4143;
	double averZ=16.4429;
	Lorentz Lor_Angle;

};


class bmt_strip
{
public:
  
  double sigma; // Transverse diffusion value computed from SigmaDrift
  double sigma_phi; // sigma/radius of the tile... for Z-detector 
  
  vector<double> FindStrip( int layer, int sector, G4ThreeVector xyz, double Edep, bmtConstants bmtc);   // Strip Finding Routine
    
  double getSigma( int layer, double x, double y, bmtConstants bmtc);     // sigma for C-detector
  int getClosestStrip( int layer, int sector, double angle, double z,bmtConstants bmtc);
  int getStripGroup(int layer, int strip, bmtConstants bmtc);
  double GetStripInfo(int layer, int sector, int strip, bmtConstants bmtc); 				   // the z position of a given C strip. Not used?
  int isInSector(int layer, double angle, bmtConstants bmtc);
  double Weight_td(int layer, int sector, int strip, double angle, double z, bmtConstants bmtc); //Compute the likelihood to get an electron
  double GetBinomial(double n, double p); //Compute the number of electrons collected following the likelihood from Weight_td
};

#endif










