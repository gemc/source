#ifndef MPrimaryGeneratorAction_h
#define MPrimaryGeneratorAction_h 1

// G4 headers
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"

// %%%%%%%%%%%%%
// For reading StdHep files
// %%%%%%%%%%%%%

#include "lStdHep.hh"
using namespace UTIL;


// gemc headers
#include "options.h"

// C++ headers
#include <fstream>


// %%%%%%%%%%%%%%%%
// Class definition
// %%%%%%%%%%%%%%%%
class MPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		MPrimaryGeneratorAction(goptions*);
		~MPrimaryGeneratorAction();
		
	public:
		void GeneratePrimaries(G4Event* anEvent);
		goptions *gemcOpt;
		double GEN_VERBOSITY;
  		double getBeamPol()
		{
			return beamPol;
		}
    
		// LUND information
		int nparticles;
		vector<double> lundUserDefined;   ///< user defined infos in the LUND header
	
	
	private:
		string input_gen;                 ///< Input Option
		string hd_msg;                    ///< Head Message Log
		G4ParticleTable* particleTable;   ///< Geant4 Particle Table
		
		// Primary Beam
		G4ParticleDefinition *Particle;   ///< Particle type
		double mom,   dmom,    Mom;       ///< beam momentum, delta momentum, randomized momentum
		double theta, dtheta,  Theta;     ///< theta, delta theta, randomized theta
		double phi,   dphi,    Phi;       ///< phi, delta phi, randomized phi
		double vx, vy, vz;                ///< Beam Vertex coordinates
		double dvr, dvz;                  ///< Deltas Beam Vertex: Radius and z-vertex
		double Vx, Vy, Vz;                ///< Randomized Beam Vertex coordinates
		double polDeg, polTheta, polPhi;  ///< Polarization degree and  direction
		G4ThreeVector beam_dir;           ///< beam direction
		G4ThreeVector beam_vrt;           ///< beam vertex
		G4ThreeVector beam_pol;           ///< beam polarization vector
		double ctheta;                    ///< customized theta direction for the z axis
		double cphi;                      ///< customize phi direction for the z axis
		
		
		// Generators Input Files
		ifstream  gif;                    ///< Generator Input File
		string    gformat;                ///< Generator Format. Supported: LUND.
		string    gfilename;              ///< Input Filename
    	double    beamPol;                ///< Beam Polarization as from the LUND format, it
	
        lStdHep   *stdhep_reader;         /// Handle to the object for reading StdHep files.
                                          ///< is only along the z axis
		
		// Luminosity Beam
		G4ParticleDefinition *L_Particle;  ///< Luminosity Particle type
		double L_mom,  L_dmom,  L_Mom;     ///< Luminosity beam momentum, delta momentum, randomized momentum
		double L_theta, L_dtheta, L_Theta; ///< Luminosity theta,  delta theta, randomized theta
		double L_phi, L_dphi, L_Phi;       ///< Luminosity phi, delta phi, randomized phi
		double L_vx, L_vy, L_vz;           ///< Luminosity Beam Vertex coordinates
		double L_dvr,  L_dvz;              ///< Luminosity Deltas Beam Vertex: Radius and z-vertex
		int NP;                            ///< Number of Luminosity Particles per event
		double TWINDOW;                    ///< Time Window
		double TBUNCH;                     ///< Time Between Bunches
		G4ThreeVector L_beam_dir;          ///< Luminosity beam direction
		G4ThreeVector L_beam_vrt;          ///< Luminosity beam vertex
		
		// Luminosity Beam2
		G4ParticleDefinition *L2_Particle;    ///< Luminosity Particle type
		double L2_mom, L2_dmom,  L2_Mom;      ///< Luminosity beam momentum, delta momentum, randomized momentum
		double L2_theta, L2_dtheta, L2_Theta; ///< Luminosity theta,  delta theta, randomized theta
		double L2_phi, L2_dphi, L2_Phi;       ///< Luminosity phi, delta phi, randomized phi
		double L2_vx, L2_vy, L2_vz;           ///< Luminosity Beam Vertex coordinates
		double L2_dvr,  L2_dvz;               ///< Luminosity Deltas Beam Vertex: Radius and z-vertex
		int NP2;                              ///< Number of Luminosity Particles per event
		double TBUNCH2;                       ///< Time Between Bunches
		G4ThreeVector L2_beam_dir;            ///< Luminosity beam direction
		G4ThreeVector L2_beam_vrt;            ///< Luminosity beam vertex
		
		G4ParticleGun* particleGun;
		void setBeam();
    
};

#endif


