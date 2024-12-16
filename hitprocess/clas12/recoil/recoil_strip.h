#ifndef recoil_strip_H
#define recoil_strip_H 1

// geant4
#include "G4ThreeVector.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;
#include <string>
#include <vector>
#include <map>

using namespace std;

//#include "Lorentz.h"

class recoilConstants
{
public:
	
	// database
	int    runNo;
	string date;
	string connection;
	char   database[80];	
	
	// trapezoid dimensions
	double Xhalf;   // half zidth
	double Yhalf;           // half height
	double Zhalf;           // Z half length
	
	// strip geometrical characteristics;
	
	int number_strip_chamber_u[3];
        int number_strip_chamber_v[3];
        int number_strip_chamber[3];
  
	double stripU_stereo_angle ; // angle between strip and trapezoid base in degree
	double stripU_pitch ;  //mm
	double stripU_width[4] ;  // mm

	double stripV_stereo_angle ; // angle between strip and trapezoid base in degree
	double stripV_pitch ;  //mm
	double stripV_width[4] ;  // mm
	
	// Detector response characteristics
	double w_i; //ionization potential assumed to be 25 eV
	double sigma_td;         // effective value to take into account transverse diffusion + charge dispersion
	int nb_sigma;            // Number of sigma to study around the closest strip
	double gain;
	
	// drift velocity
	double v_drift; // velocity drift [cm/ns]
	double sigma_time ; // time resolution 20 ns
	
	double v_eff_readout = 0; // effective signal velocity
private:
	double stereo_angle;
	double strip_pitch;
	double strip_width[4];
	string kind_strip;
	
public:

	inline void get_strip_info(string a){
		if (a == "strip_u") {
			stereo_angle = stripU_stereo_angle;
			strip_pitch = stripU_pitch;
			copy(std::begin(stripU_width), std::end(stripU_width), std::begin(strip_width));
			kind_strip = a;
		}
		if (a == "strip_v") {
			stereo_angle = stripV_stereo_angle;
			strip_pitch = stripV_pitch;
			copy(std::begin(stripV_width), std::end(stripV_width), std::begin(strip_width));			
			kind_strip = a;
		}
	}

	inline double get_stereo_angle() {return stereo_angle;}
	inline double get_strip_pitch() {return strip_pitch;}
	inline double get_strip_width(int strip) {
		double width;
		width = strip_width[0];
		
		return width;
	}
	inline string get_strip_kind() {return kind_strip;}
};

class recoil_strip_found
{
public:
	int numberID;
	double weight;
	double time;
};

class recoil_strip
{
public:
	
  int number_of_strip;
	int Nel;                 // number of electrons (Nt) for a given hit
	double y_real; 		 // y position of the hit after transverse diffusion
	double x_real;		 // x position of the hit after transverse diffusion
	double strip_x;          // strip_x is the middle of the strips
	double strip_y;          // strip_y is the position of the strips
	double strip_length;     // length of the strip
	G4ThreeVector strip_endpoint1;
	G4ThreeVector strip_endpoint2;
	
	vector<recoil_strip_found> FindStrip( G4ThreeVector   xyz , double Edep, recoilConstants recoilc, double time); // Strip Finding Routine
	void Carac_strip(int strip, recoilConstants recoilc); //length of the strip
	double Weight_td(int strip, double x, double y, double z, recoilConstants recoilc); //Compute the fraction of Nel falling onto the strip, depending on x,y in the FMT coordinate system
	bool Build_strip(int strip, recoilConstants recoilc );
	
	G4ThreeVector intersectionPoint(double m, double c, G4ThreeVector A, G4ThreeVector B);
	double cal_length(G4ThreeVector A, G4ThreeVector B);
	bool pointOnsegment(G4ThreeVector X, G4ThreeVector A, G4ThreeVector B);
	G4ThreeVector change_of_coordinates(G4ThreeVector   xyz , recoilConstants recoilc);
	
	double GetBinomial(double n, double p);
        int Number_of_strip(recoilConstants recoilc);
	int strip_id(int i, recoilConstants recoilc);
	
};



#endif
