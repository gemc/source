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

	double SigmaDrift;	 // Max transverse diffusion value (GEMC value)
	double hDrift;			 // Size of the drift gap
	double hStrip2Det;	 // Distance between strips and the middle of the conversion gap (~half the drift gap)
	double ThetaL;		    // Lorentz angle
	double fieldScale; 	 // Scaling of the field - 1 means a 5 Tesla field along z

	// THE GEOMETRY CONSTANTS
	const static int NREGIONS = 3  ;	// 3 regions of MM

	// Z detector characteristics
	double CRZRADIUS[NREGIONS] ; 	       // the radius of the Z detector in mm
	int CRZNSTRIPS[NREGIONS] ; 	       // the number of strips
	double CRZSPACING[NREGIONS] ;        // the strip spacing in mm
	double CRZWIDTH[NREGIONS] ; 	       // the strip width in mm
	double CRZLENGTH[NREGIONS] ; 	       // the strip length in mm
	double CRZZMIN[NREGIONS] ; 	       // PCB upstream extremity mm
	double CRZZMAX[NREGIONS] ; 	       // PCB downstream extremity mm
	double CRZOFFSET[NREGIONS] ; 	       // Beginning of strips in mm
	double CRZXPOS[NREGIONS]; 	          // Distance on the PCB between the PCB first edge and the edge of the first strip in mm
	double CRZEDGE1[NREGIONS][NREGIONS]; // the angle of the first edge of each PCB detector A, B, C
	double CRZEDGE2[NREGIONS][NREGIONS]; // the angle of the second edge of each PCB detector A, B, C

	// C detector characteristics
	double CRCRADIUS[NREGIONS]; 		    // the radius of the Z detector in mm
	int CRCNSTRIPS[NREGIONS]; 			    // the number of strips
	double CRCSPACING[NREGIONS]; 		    // the strip spacing in mm
	double CRCLENGTH[NREGIONS]; 		    // the strip length in mm
	double CRCZMIN[NREGIONS]; 			    // PCB upstream extremity mm
	double CRCZMAX[NREGIONS]; 			    // PCB downstream extremity mm
	double CRCOFFSET[NREGIONS]; 		    // Beginning of strips in mm
	double CRCXPOS[NREGIONS]; 			    // Distance on the PCB between the PCB first edge and the edge of the first strip in mm
	double CRCEDGE1[NREGIONS][NREGIONS]; // the angle of the first edge of each PCB detector A, B, C
	double CRCEDGE2[NREGIONS][NREGIONS]; // the angle of the second edge of each PCB detector A, B, C
	vector<vector<int> >     CRCGROUP; 	 // Number of strips with same width
	vector<vector<double> >  CRCWIDTH;   // the width of the corresponding group of strips

	void changeFieldScale(double newFieldScale)
	{
		fieldScale = newFieldScale;
		ThetaL     = fieldScale*20.*degree;
	}

};


class bmt_strip
{
public:



	vector<double> FindStrip( int layer, int sector, G4ThreeVector xyz, double Edep, bmtConstants bmtc);   // Strip Finding Routine


	double getEnergyFraction(double z0, double z, double sigma); 			// gaussian pdf
	int getDetectorIndex(int sector);										// index of the detector (0...2): sector 1 corresponds to detector B, 2 to A,  3 to C in the geometry created by GEMC

	double getSigmaLongit( int layer, double x, double y, bmtConstants bmtc);     // sigma for C-detector
	double getSigmaAzimuth(int layer, double x, double y, bmtConstants bmtc);     // sigma for Z-detectors
	int getZStrip(int layer, double angle, bmtConstants bmtc); 						   // the Z strip as a function of azimuthal angle
	int getCStrip(int layer, double trk_z, bmtConstants bmtc); 					 	   // the Z strip as a function of z
	double CRCStrip_GetZ(int layer, int strip, bmtConstants bmtc); 				   // the z position of a given C strip. Not used?
	double CRZStrip_GetPhi(int sector, int layer, int strip, bmtConstants bmtc);	// the phi angle of a given Z strip. Not used?
	int isInSector(int layer, double angle, bmtConstants bmtc);							// the sector according to this geometry defined in fillinfos()


};

#endif










