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

public:
	~hipo_output(){;}  ///< event is deleted in WriteEvent routine
	static outputFactory *createOutput() {return new hipo_output;}
	
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
	
};


#endif
