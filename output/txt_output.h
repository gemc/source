/// \file txt_output.h 
/// Defines the Text Output Class.\n
/// The pointer to the output stream
/// is passed by the outputContainer class.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n
#ifndef TXT_OUTPUT_H
#define TXT_OUTPUT_H 1



// gemc headers
#include "outputFactory.h"

// Class definition
class txt_output : public outputFactory
{
	public:
		~txt_output(){;}  ///< event is deleted in WriteEvent routine
		static outputFactory *createOutput() {return new txt_output;}

	// record the simulation conditions on the file
	void recordSimConditions(outputContainer*, map<string, string>);

	// write header
	void writeHeader(outputContainer*, map<string, double>, gBank);

	// write generated particles
	void writeGenerated(outputContainer*, vector<generatedParticle>, map<string, gBank> *banksMap);
	
	// format output and set insideBank
	void initBank(outputContainer*, gBank);
	
	// write geant4 raw integrated info
	void writeG4RawIntegrated(outputContainer*, vector<hitOutput>,  string, map<string, gBank>*);
		
	// write geant4 digitized integrated info
	void writeG4DgtIntegrated(outputContainer*, vector<hitOutput>,  string, map<string, gBank>*);

	// write geant4 true info for every step
	virtual void writeG4RawAll(outputContainer*, vector<hitOutput>, string, map<string, gBank>*);


	// write event and close stream if necessary
	// nothing to be done for txt
	void writeEvent(outputContainer*);


	map<string, bool> insideBank;
};
#endif
