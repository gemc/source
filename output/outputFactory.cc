// gemc headers
#include "outputFactory.h"
#include "options.h"
#include "string_utilities.h"
#include "evio_output.h"
#include "txt_output.h"

// C++ headers
#include <fstream>
using namespace std;

#define EVIO_BUFFER 30000000

outputFactory *getOutputFactory(map<string, outputFactoryInMap> *outputFactoryMap, string outputType)
{
	
	if(outputFactoryMap->find(outputType) == outputFactoryMap->end())
	{
		cout << endl << endl << "  >>> WARNING: Output type <" << outputType << "> NOT FOUND IN  Output Map." << endl;
		return NULL;
	}
	
	return (*outputFactoryMap)[outputType]();
}

outputContainer::outputContainer(goptions Opts)
{
	// EVIO Buffer size set to 30M words
	int evio_buffer = EVIO_BUFFER;
	
	gemcOpt = Opts;
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Output File: >> ";
	
	string optf = gemcOpt.optMap["OUTPUT"].args;
	outType.assign(optf, 0, optf.find(",")) ;
	outFile.assign(optf,    optf.find(",") + 1, optf.size()) ;
	
	if(outType != "no") cout << hd_msg << " Opening output file \"" << TrimSpaces(outFile) << "\"." << endl;
	if(outType == "txt")  txtoutput = new ofstream(TrimSpaces(outFile).c_str());
	if(outType == "evio")
	{
		pchan = new evioFileChannel(TrimSpaces(outFile).c_str(), "w", evio_buffer);
		pchan->open();
	}
}

outputContainer::~outputContainer()
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Output File: >> ";
	
	if(outType != "no")   cout << " Closing " << outFile << "." << endl;
	if(outType == "txt")  txtoutput->close();
	if(outType == "evio")
	{
		pchan->close();
		delete pchan;
	}
}


map<string, outputFactoryInMap> registerOutputFactories()
{
	map<string, outputFactoryInMap> outputMap;
	
	outputMap["txt"]   =   &txt_output::createOutput;
	outputMap["evio"]  =  &evio_output::createOutput;
	
	return outputMap;
}

int generatedParticle::getVariableFromStringI(string what)
{
	if(what == "pid") return PID;
	return 0;
}

double generatedParticle::getVariableFromStringD(string what)
{
	if(what == "px") return momentum.x();
	if(what == "py") return momentum.y();
	if(what == "pz") return momentum.z();
	if(what == "vx") return vertex.x();
	if(what == "vy") return vertex.y();
	if(what == "vz") return vertex.z();

	return 0;
}
