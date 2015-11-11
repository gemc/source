#ifndef bmt_strip_H
#define bmt_strip_H 1


#include <vector>
using namespace std;

/// \class bmt_strip
/// <b> bmt_strip </b>\n\n
/// Micromegas strip finding routine\n


class bmt_strip
{
	public:

		// THE GEOMETRY CONSTANTS
		const static int NREGIONS = 3  ;	// 3 regions of MM
		
		//Z detector characteristics
		double CRZRADIUS[NREGIONS] ; 		// the radius of the Z detector in mm
		int CRZNSTRIPS[NREGIONS] ; 			// the number of strips
		double CRZSPACING[NREGIONS] ; 		// the strip spacing in mm
		double CRZWIDTH[NREGIONS] ; 		// the strip width in mm
		double CRZLENGTH[NREGIONS] ; 		// the strip length in mm
		double CRZZMIN[NREGIONS] ; 			// PCB upstream extremity mm
		double CRZZMAX[NREGIONS] ; 			// PCB downstream extremity mm
		double CRZOFFSET[NREGIONS] ; 		// Beginning of strips in mm
		vector<vector<double> >  CRZEDGE1; 	// the angle of the first edge of each PCB detector A, B, C
		vector<vector<double> >  CRZEDGE2; 	// the angle of the second edge of each PCB detector A, B, C
		double CRZXPOS[NREGIONS]; 			// Distance on the PCB between the PCB first edge and the edge of the first strip in mm

		//C detector characteristics
		double CRCRADIUS[NREGIONS]; 		// the radius of the Z detector in mm
		int CRCNSTRIPS[NREGIONS]; 			// the number of strips
		double CRCSPACING[NREGIONS]; 		// the strip spacing in mm
		double CRCLENGTH[NREGIONS]; 		// the strip length in mm
		double CRCZMIN[NREGIONS]; 			// PCB upstream extremity mm
		double CRCZMAX[NREGIONS]; 			// PCB downstream extremity mm
		double CRCOFFSET[NREGIONS]; 		// Beginning of strips in mm
		vector<vector<int> >  CRCGROUP; 	// Number of strips with same width
		vector<vector<double> >  CRCWIDTH;	// the width of the corresponding group of strips
		vector<vector<double> >  CRCEDGE1; 	// the angle of the first edge of each PCB detector A, B, C
		vector<vector<double> >  CRCEDGE2; 	// the angle of the second edge of each PCB detector A, B, C
		double CRCXPOS[NREGIONS]; 			// Distance on the PCB between the PCB first edge and the edge of the first strip in mm

		// THE SIMULATION CONSTANTS
		double Pi ;
		double SigmaDrift;
		double hDrift;
		double hStrip2Det;
		double ThetaL ;

		void fill_infos(); 
		vector<double> FindStrip( int layer, int sector, double x, double y, double z, double Edep);   // Strip Finding Routine
		double toRadians(double angleDegrees);
		bool isInFiducial(int sector, int layer, double[] x);  // a boolean to indicate if the hit is within sensitive area
		double getSigmaAzimuth(int layer, double x, double y); // sigma for Z-detectors
		double getSigmaLongit(int layer, double x, double y);  // sigma for C-detector
		vector<double> smearedPosition(int layer, double x, double y, double z); //diffused position
		int getZStrip(int layer, double angle); //the Z strip as a function of azimuthal angle
		double CRCStrip_GetZ(int sector, int layer, int strip); // the z position of a given C strip
};

#endif
