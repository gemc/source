/// \file hipo_output.h
/// Defines the Hipo Output Class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef HIPO_OUTPUT_H
#define HIPO_OUTPUT_H 1

/*

	Hipo Library: https://github.com/gavalian/hipo
	CLAS12 Bank Definitions: https://github.com/JeffersonLab/clas12-offline-software/tree/development/etc/bankdefs/hipo4

	The Digitized banks are in data.json

 */


// gemc headers
#include "outputFactory.h"


// Class definition
class hipo_output : public outputFactory
{
private:
	static map<string, double> fieldScales;
	static int rasterInitialized;
	static double rasterP0[2];
	static double rasterP1[2];

public:
	~hipo_output(){
		if(outEvent) {
			delete outEvent;
		}
		if(trueInfoBank) {
			delete trueInfoBank;
		}

	}  ///< event is deleted in WriteEvent routine
	static outputFactory *createOutput() {return new hipo_output;}
	
	// prepare event
	virtual void prepareEvent(outputContainer* output, map<string, double> *configuration) ;

	// record the simulation conditions on the file
	void recordSimConditions(outputContainer*, map<string, string>);
	
	// write header
	void writeHeader(outputContainer*, map<string, double>, gBank);
	
	// write user infos header
	void writeUserInfoseHeader(outputContainer*, map<string, double>);
	
	// format output and set insideBank
	void initBank(outputContainer*, gBank, int what);
	
	// write generated particles
	void writeGenerated(outputContainer*, vector<generatedParticle>, map<string, gBank> *banksMap, vector<userInforForParticle> userInfo);
	
	// write ancestors
	virtual void writeAncestors (outputContainer*, vector<ancestorInfo>, gBank);
	
	// write RF Signal
	virtual void writeRFSignal(outputContainer*, FrequencySyncSignal, gBank);
	
	// write geant4 raw integrated info
	void writeG4RawIntegrated(outputContainer*, vector<hitOutput>,  string, map<string, gBank>*);
	
	// write geant4 digitized integrated info
	void writeG4DgtIntegrated(outputContainer*, vector<hitOutput>,  string, map<string, gBank>*);
	
	// write geant4 charge / time (as seen by electronic) info
	virtual void writeChargeTime(outputContainer*, vector<hitOutput>, string, map<string, gBank>*);
	
	// write geant4 true info for every step
	virtual void writeG4RawAll(outputContainer*, vector<hitOutput>, string, map<string, gBank>*);
	
	// write fadc mode 1 (full signal shape) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
	virtual void writeFADCMode1(outputContainer*, vector<hitOutput>, int);
	
	// write fadc mode 1 (full signal shape) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
	// This method should be called once at the end of event action, and the 1st argument 
	// is a map<int crate_id, vector<hitoutput> (vector of all hits from that crate) >
	virtual void writeFADCMode1( map<int, vector<hitOutput> >, int);
	
	// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
	virtual void writeFADCMode7(outputContainer*, vector<hitOutput>, int);
	
	// write event and close stream if necessary
	void writeEvent(outputContainer*) ;

	hipo::event *outEvent = nullptr;
	hipo::bank *trueInfoBank = nullptr;

	// needed to correctly index
	int lastHipoTrueInfoBankIndex  = 0;

	// map from hittype to hipo detector id
	// notice the hittype is gemc specific
	// need to change naming convention here to match reconstruction?
	// defined here: https://github.com/JeffersonLab/clas12-offline-software/blob/8ed53986f8b1a2e6f3c5a63b1e6f6d7fd88020c9/common-tools/clas-detector/src/main/java/org/jlab/detector/base/DetectorType.java
	map<string, int> detectorID = {
		{"bmt",     1},
		{"bst",     2},
		{"cnd",     3},
		{"ctof",    4},
		{"dc",      6},
		{"ecal",    7},
		{"fmt",     8},
		{"ft_cal",  10},
		{"ft_hodo", 11},
		{"ft_trk",  13},
		{"ftof",    12},
		{"htcc",    15},
		{"ltcc",    16},
		{"rich",    18},
		{"rtpc",    19},
		{"band",    21},
		{"urwell",  23},
		{"atof",    24},
		{"ahdc",    25},
		{"flux",   100}
	};

	// returns detectorID from map, given hitType
	int getDetectorID(string hitType) ;

	// true info variable names are changed a bit in hipo schema
	map<string, string> trueInfoNamesMap = {
		{"pid",     "pid"},
		{"mpid",    "mpid"},
		{"tid",     "tid"},
		{"mtid",    "mtid"},
		{"otid",    "otid"},
		{"trackE",  "trackE"},
		{"totEdep", "totEdep"},
		{"avg_x",   "avgX"},
		{"avg_y",   "avgY"},
		{"avg_z",   "avgZ"},
		{"avg_lx",  "avgLx"},
		{"avg_ly",  "avgLy"},
		{"avg_lz",  "avgLz"},
		{"px",      "px"},
		{"py",      "py"},
		{"pz",      "pz"},
		{"vx",      "vx"},
		{"vy",      "vy"},
		{"vz",      "vz"},
		{"mvx",     "mvx"},
		{"mvy",     "mvy"},
		{"mvz",     "mvz"},
		{"avg_t",   "avgT"},
		{"nsteps",  "nsteps"},
		{"procID",  "procID"},
		{"hitn",    "hitn"}

	};

	// returns hipo name from true info var name
	string getHipoVariableName(string trueInfoVar) ;
	
	
	

};


#endif
