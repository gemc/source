/// \file run_conditions.h
/// Defines the Run Conditions class.\n
/// The Run Conditions are defined in an XML file called gcard\n
/// Run condition includes:
/// - detectory system and their factory\n
/// - tilts and displacements
/// - particle generation definitions
/// - options defined in the gemc option map
/// The gemc option map is filled with the run conditions
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef run_conditions_H
#define run_conditions_H

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

// Qt headers
#include <QDomDocument>

// gemc headers
#include "options.h"

// C++ headers
#include <iostream>
using namespace std;


/// \class detectorCondition
/// <b> detectorCondition </b>\n\n
/// This is the detector Conditions class.\n
/// For each detector its shifts and rotation relative to the nominal position are recorded here. \n
/// The factory type of the detector is also recorded here.
class detectorCondition
{
	private:
		G4ThreeVector    pos;  // Shift relative to the nominal position
		G4RotationMatrix rot;  // Rotation relative to the nominal position
		G4ThreeVector    vrot; // Rotation Vector (ordered X,Y,Z)

		bool presentFlag;      // by default this is false. if existance is set in the gcard, this is set to true
		int is_present;        // by default set to 1. Can be set to 0 to remove a detector from the simulation
		string system;         // detector system
		string factory;        // factory that builds the detector
		string variation;      // variation of the detector. Default is "main"
		int    run_number;     // Run Number selected for this detector
		
	public:
	detectorCondition(){is_present = 0; presentFlag = false;}
		detectorCondition(string f)
		{
			factory    = f;
			is_present = 0;
			variation  = "main";
			run_number = 1;
			presentFlag = false;
		}
		~detectorCondition(){;}
		
		void set_position(string X, string Y, string Z);
		void set_rotation(string X, string Y, string Z);
		void set_existance(string exist);
		void set_factory(string f)   {factory = f;}
		void set_variation(string v) {variation = v;}
		void set_system(string s)    {system = s;}
		void set_run_number(int r)   {run_number = r;}
		
		G4ThreeVector    get_position()  {return pos;}
		G4ThreeVector    get_vrotation() {return vrot;}
		G4RotationMatrix get_rotation()  {return rot;}
		int              get_existance() {if(presentFlag) return is_present; else return 2;} // 2 means existance is not set
		string           get_factory()   {return factory;}
		string           get_variation() {return variation;}
		string           get_system()    {return system;}
		int              get_run_number(){return run_number;}
};


/// \class runConditions
/// <b> runConditions </b>\n\n
/// This is the gemc Run Conditions class. It contains a map of detectorCondition
/// and all options from the gcard file.\n
class runConditions
{
	public:
		runConditions(goptions);
		runConditions(){;}
		~runConditions();
		
		// Map of detectorCondition. Map Key = detector name.
		map<string, detectorCondition> detectorConditionsMap;
		
		// Returns a map<string, string> with the detector systems present in gemc
		// This map is then written in the output stream
		map<string, string> getDetectorConditionsMap();
		
		
		int get_run_number(string detector);
		string get_variation(string detector);
		string get_system(string detector);
		
		map<string, string> get_systems();
	
	private:
		QDomDocument domDocument;
		
};

int check_if_factory_is_needed(map<string, detectorCondition>, string);



class runWeights
{
	public:
		runWeights(goptions);
		runWeights(){;}
    	~runWeights(){;}

		int runNo;
		int getRunNumber(int n);
		bool isNewRun;
		int defaultRunNumber;
	
	private:
		// map with weights as coming from the file
		map<int, double> w;
		
		// map with numnber of events for each run, based on weight map
		map<int, int> n;
	
		int startEvent;
	
	
};




#endif






