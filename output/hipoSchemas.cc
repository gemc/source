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

	// CLAS12 schemas are defined at https://github.com/JeffersonLab/clas12-offline-software/tree/development/etc/bankdefs/hipo4
	// Detectors: https://github.com/JeffersonLab/clas12-offline-software/blob/development/etc/bankdefs/hipo4/data.json
	runConfigSchema = hipo::schema("RUN::config", 10000, 11);

	// detectors
	bmtADCSchema    = hipo::schema("BMT::adc",    20100, 11);
	bstADCSchema    = hipo::schema("BST::adc",    20200, 11);
	cndADCSchema    = hipo::schema("CND::adc",    20300, 11);
	cndTDCSchema    = hipo::schema("CND::tdc",    20300, 12);
	ctofADCSchema   = hipo::schema("CTOF::adc",   20400, 11);
	ctofTDCSchema   = hipo::schema("CTOF::tdc",   20400, 12);
	dcTDCSchema     = hipo::schema("DC::tdc",     20600, 12);
	dcDOCASchema    = hipo::schema("DC::doca",    20600, 14);
	ecalADCSchema   = hipo::schema("ECAL::adc",   20700, 11);
	ecalTDCSchema   = hipo::schema("ECAL::tdc",   20700, 12);
	fmtADCSchema    = hipo::schema("FMT::adc",    20800, 11);
	ftcalADCSchema  = hipo::schema("FTCAL::adc",  21000, 11);
	fthodoADCSchema = hipo::schema("FTHODO::adc", 21100, 11);
	ftofADCSchema   = hipo::schema("FTOF::adc",   21200, 11);
	ftofTDCSchema   = hipo::schema("FTOF::tdc",   21200, 12);
	ftrkTDCSchema   = hipo::schema("FTTRK::adc",  21300, 11);
	htccADCSchema   = hipo::schema("HTCC::adc",   21500, 11);
	htccTDCSchema   = hipo::schema("HTCC::tdc",   21500, 12);
	ltccADCSchema   = hipo::schema("LTCC::adc",   21600, 11);
	ltccTDCSchema   = hipo::schema("LTCC::tdc",   21600, 12);
	rfADCSchema     = hipo::schema("RF::adc",     21700, 11);
	rfTDCSchema     = hipo::schema("RF::adc",     21700, 12);
	richTDCSchema   = hipo::schema("RICH::tdc",   21800, 12);
	rtpcADCSchema   = hipo::schema("RTPC::adc",   21900, 11);
	rtpcPOSSchema   = hipo::schema("RTPC::pos",   21900, 14);
	bandADCSchema   = hipo::schema("BAND::adc",   22100, 11);
	bandTDCSchema   = hipo::schema("BAND::tdc",   22100, 12);
	helADCSchema    = hipo::schema("HEL::adc",    22000, 11);
	helFLIPSchema   = hipo::schema("HEL::flip",   22000, 12);
	helONLINESchema = hipo::schema("HEL::online", 22000, 13);
	rawADCSchema    = hipo::schema("RAW::adc",    20000, 11);
	rawTDCSchema    = hipo::schema("RAW::tdc",    20000, 12);
	rawSCALERSchema = hipo::schema("RAW::scaler", 20000, 13);
	rawVTPSchema    = hipo::schema("RAW::vtp",    20000, 14);
	rawEPICSSchema  = hipo::schema("RAW::epics",  20000, 15);


	// Defining structure of the schema (bank)
	// The columns in the banks (or leafs, if you like ROOT)
	// are given comma separated with type after the name.
	// Available types : I-integer, S-short, B-byte, F-float,
	//                   D-double, L-long

	runConfigSchema.parse("run/I, event/I, unixtime/I, trigger/L, timestamp/L, type/B,mode/B, torus/F, solenoid/F");

	// detectors
	bmtADCSchema.parse(    "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");
	bstADCSchema.parse(    "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, timestamp/L");
	ctofADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	ctofTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	fmtADCSchema.parse(    "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");

	// todolater
	cndADCSchema.parse(    "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	cndTDCSchema.parse(    "sector/B, layer/B, component/S, order/B, TDC/I");

	dcTDCSchema.parse(     "sector/B, layer/B, component/S, order/B, TDC/I");

	ecalADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	ecalTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	ftcalADCSchema.parse(  "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	fthodoADCSchema.parse( "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	ftofADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	ftofTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	ftrkTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S, integral/I, timestamp/L");
	htccADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	htccTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	ltccADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	ltccTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	rfADCSchema.parse(     "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	rfTDCSchema.parse(     "sector/B, layer/B, component/S, order/B, TDC/I");
	richTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	rtpcADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	rtpcPOSSchema.parse(   "step/I, time/F, energy/F, posx/F, posy/F, posz/F, phi/F, tid/F");
	bandADCSchema.parse(   "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	bandTDCSchema.parse(   "sector/B, layer/B, component/S, order/B, TDC/I");
	helADCSchema.parse(    "sector/B, layer/B, component/S, order/B, ADC/I, time/F, ped/S");
	helFLIPSchema.parse(   "run/I, event/I, timestamp/L, helicity/B, helicityRaw/B, pair/B, pattern/B, status/B");
	helONLINESchema.parse( "helicity/B, helicityRaw/B");
	rawADCSchema.parse(    "crate/B, slot/B, channel/S, order/B, ADC/I, time/F, ped/S");
	rawTDCSchema.parse(    "crate/B, slot/B, channel/S, order/B, TDC/I");
	rawSCALERSchema.parse( "crate/B, slot/B, channel/S, helicity/B, quartet/B, value/L");
	rawVTPSchema.parse(    "crate/B, word/I");
	rawEPICSSchema.parse(  "json/B");
	emptySchema.parse(     "empty/B");

	schemasToLoad["RUN::config"] = runConfigSchema;

	// Central Detector
	schemasToLoad["BMT::adc"]    = bmtADCSchema;
	schemasToLoad["BST::adc"]    = bstADCSchema;
	schemasToLoad["CTOF::adc"]   = ctofADCSchema;
	schemasToLoad["CTOF::tdc"]   = ctofTDCSchema;
	schemasToLoad["FMT::adc"]    = fmtADCSchema;
	schemasToLoad["DC::tdc"]     = dcTDCSchema;


	cout << " Done defining Hipo4 schemas." << endl;

}

#include <cstring>
// type: 0 = adc, 1 = tdc
hipo::schema HipoSchema :: getSchema(string schemaName, int type) {

	string schemaType = type == 0 ? "adc" : "tdc";

	string toUpperS = schemaName;
	transform(toUpperS.begin(), toUpperS.end(), toUpperS.begin(), ::toupper);
	string thisSchema = toUpperS + "::" + schemaType;

	if(schemasToLoad.find(thisSchema) != schemasToLoad.end() ) {
		return schemasToLoad[thisSchema];
	} else {
		return emptySchema;
	}
}


void outputContainer::initializeHipo(string outputFile) {

	hipoSchema = new HipoSchema();

	cout << " Initializing hipo4 writer to filename:" << outputFile << endl;
	hipoWriter = new hipo::writer();

	// Open a writer and register schemas with the writer.
	// The schemas have to be added to the writer before openning
	// the file, since they are written into the header of the file.
	for (auto &schema: hipoSchema->schemasToLoad) {
		hipoWriter->getDictionary().addSchema(schema.second);
	}
	
	hipoWriter->open(outputFile.c_str());

}
