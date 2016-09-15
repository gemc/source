/// \file evio_output.h 
/// Defines the EVIO Output Class.\n
/// The pointer to the evioFileChannel
/// is passed by the outputContainer class.\n
/// The evioDOMTree and evioDOMNodeP are 
/// created on the heap.
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef EVIO_OUTPUT_H
#define EVIO_OUTPUT_H 1

// gemc headers
#include "outputFactory.h"


// Class definition
class evio_output : public outputFactory
{
	public:
		~evio_output(){;}  ///< event is deleted in WriteEvent routine
		static outputFactory *createOutput() {return new evio_output;}
  
	// record the simulation conditions on the file
	void recordSimConditions(outputContainer*, map<string, string>);

	// write header
	void writeHeader(outputContainer*, map<string, double>, gBank);

	// format output and set insideBank
	void initBank(outputContainer*, gBank, int what);
 
	// write generated particles
	void writeGenerated(outputContainer*, vector<generatedParticle>, map<string, gBank> *banksMap);

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
	virtual void writeFADCMode1(outputContainer*, vector<hitOutput>);

	// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
	virtual void writeFADCMode7(outputContainer*, vector<hitOutput>);

	// write event and close stream if necessary
	void writeEvent(outputContainer*) ;

	
	evioDOMTree *event;
	evioDOMNodeP detectorBank;
	map<string, evioDOMNodeP> detectorRawIntBank; 
	map<string, evioDOMNodeP> detectorDgtIntBank;
	map<string, evioDOMNodeP> detectorRawAllBank;
	map<string, evioDOMNodeP> detectorChargeTimeBank;

	map<string, bool> insideBank;
	map<string, bool> insideRawIntBank;
	map<string, bool> insideDgtIntBank;
	map<string, bool> insideRawAllBank;
	map<string, bool> insideChargeTimeBank;

	int evn;
	
};

// returns a evioDOMNodeP based on the type specified by the string
// only double is casted back
evioDOMNodeP addVariable(int, int, string, double);
evioDOMNodeP addVariable(int, int, string, int);
evioDOMNodeP addVariable(int, int, string, string);

evioDOMNodeP addVector(int, int, string, vector<double>);
evioDOMNodeP addVector(int, int, string, vector<int>);
evioDOMNodeP addVector(int, int, string, vector<string>);


#endif
