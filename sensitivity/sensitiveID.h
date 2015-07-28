/// \file sensitiveID.h
/// Defines the gemc Sensitive Identifier class\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef sensitiveID_H
#define sensitiveID_H 1

// gemc headers
#include "options.h"

// C++ headers
#include <string>
#include <vector>
using namespace std;


/// \class sensitiveID
/// <b> sensitiveID</b>\n\n
/// The sensitiveID class contains:\n
/// - the bank tag (unique id that identifies the detector)
/// - the identifier system
/// - the signal threshold to be recorded in the output
/// -
class sensitiveID
{
	public:
		string         name;             ///< Sensitive Detector name. This has to match the bank name
		string         description;      ///< Sensitive Detector description
		vector<string> identifiers;      ///< vector of strings that uniquely identify the detector element
		double         signalThreshold;  ///< Minimum energy of the hit to be recorded in the output stream
		double         timeWindow;       ///< If two steps happens within the same TimeWindow, they belong to the same Hit
		double         prodThreshold;    ///< Geant4 Production Threshold in the detector
		double         maxStep;          ///< Geant4 Maximum Acceptable Step in the detector
		double         riseTime;         ///< rise time of the PMT signal
		double         fallTime;         ///< fall time of the PMT signal
		double         mvToMeV;          ///< from MeV to mV constant
		double         pedestal;         ///< pedestal
		double         delay;            ///< time from PMT face to signal
		string         thisFactory;      ///< Factory used to generate the sensitive detector
		
		// class constructor
		sensitiveID(string name, goptions, string factory, string variation, string system);
		sensitiveID(){;}

		friend ostream &operator<<(ostream &stream, sensitiveID SD);       ///< Overloaded "<<" for the class 'sensitiveID'
};


#endif










