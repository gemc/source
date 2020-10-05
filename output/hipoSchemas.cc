// gemc headers
#include "outputFactory.h"
#include "hipoSchemas.h"

HipoSchema :: HipoSchema()
{

	cout << " Defining Hipo4 schemas..." << endl;
	//---------------------------------------------------------
	// Define a Schema for particle bank and detector bank
	// Schema constructor takes:
	// name - name of the schema (bank)
	// groupid - a 16 bit number which identifies a group.
	//           this way banks can be grouped in skimming
	// itemid  - and item number (8 bit) to just have a
	//           unique number associated with bank.


	runConfig = hipo::schema("RUN::config", 10000, 11);




	// Defining structure of the schema (bank)
	// The columns in the banks (or leafs, if you like ROOT)
	// are given comma separated with type after the name.
	// Available types : I-integer, S-short, B-byte, F-float,
	//                   D-double, L-long


	//	"name": "RUN::config",
	//	"group": 10000,
	//	"item" : 11,
	//	"info": "Run Configuration",
	//	"entries": [
	//					{"name":"run",          "type":"I", "info":"RUN number from CODA or GEMC"},
	//					{"name":"event",        "type":"I", "info":"Event number"},
	//					{"name":"unixtime",     "type":"I", "info":"Unix time (seconds)"},
	//					{"name":"trigger",      "type":"L", "info":"trigger bits"},
	//					{"name":"timestamp",    "type":"L", "info":"time stamp from Trigger Interface (TI) board (4 nanoseconds)"},
	//					{"name":"type",         "type":"B", "info":"type of the run"},
	//					{"name":"mode",         "type":"B", "info":"run mode"},
	//					{"name":"torus",        "type":"F", "info":"torus setting relative value(-1.0 to 1.0)"},
	//					{"name":"solenoid",     "type":"F", "info":"solenoid field setting (-1.0 to 1.0)"}
	//					]
	runConfig.parse("run/I, event/I, unixtime/I, trigger/L, timestamp/L, type/B, mode/B, torus/F, solenoid/F");

	cout << " Done defining Hipo4 schemas." << endl;

}

void outputContainer::initializeHipo(string outputFile) {


	HipoSchema *hipoSchema = new HipoSchema();

	cout << " Initializing hip4 writer to filename:" << outputFile << endl;
	writer = new hipo::writer();

	// Open a writer and register schemas with the writer.
	// The schemas have to be added to the writer before openning
	// the file, since they are written into the header of the file.
	writer->getDictionary().addSchema(hipoSchema->runConfig);

	
	writer->open(outputFile.c_str());

}
