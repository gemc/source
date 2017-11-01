// gemc headers
#include "evio_output.h"
#include "utils.h"

// C++ headers
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// evio
#include "evioBankUtil.h"

#define MAXEVIOBUF 10000000
static unsigned int buf[MAXEVIOBUF];


const int evio_output::fadc_mode1_banktag = 0xe101;


// This variable is just a flag for checking, whether the FADC configuration parameters
// are written into. This should be written only once, in the 1st physics event
// it is now set false, and when configuration parameters are written, it will be
// changed to true
//bool evio_output::is_conf_written = false;


vector<int> evio_output::detector_crates = {1, 7, 13, 19, 25, 31,     // EC
	3, 9, 15, 21, 26, 33,     // PCal
	59,                       // HTCC, CTOF
        70, 71,                   // FTCalo
        72,                       // FT Hodo
        5, 11, 17, 23, 29, 35,     // FTOF        
};

// record the simulation conditions
// the format is a string for each variable
// the num is 0 for the mother bank, 1 for the data
void evio_output :: recordSimConditions(outputContainer* output, map<string, string> sims)
{
	vector<string> data;

	// writing both key and argument as one string
	for(map<string, string>::iterator it=sims.begin(); it!=sims.end(); it++)
	{
		data.push_back(it->first + ":  " + it->second + "  ");
	}

	evioDOMTree *conditionsBank = new evioDOMTree(SIMULATION_CONDITIONS_BANK_TAG, 0);
	*conditionsBank << evioDOMNode::createEvioDOMNode(SIMULATION_CONDITIONS_BANK_TAG, 1, data);

	output->pchan->write(*conditionsBank);
	delete conditionsBank;
}

// instantiates the DOM tree
// write header bank
// each variable is a domnode
void evio_output :: writeHeader(outputContainer* output, map<string, double> data, gBank bank)
{
	event = new evioDOMTree(1, 0);

	evioDOMNodeP headerBank = evioDOMNode::createEvioDOMNode(HEADER_BANK_TAG, 0);

	// timestamp
	string time = timeStamp();
	*headerBank << addVariable(HEADER_BANK_TAG, bank.getVarId("time"), "s", time);

	for(map<string, double> :: iterator it = data.begin(); it != data.end(); it++) {

		int bankId = bank.getVarId(it->first);

		// bankID 1 is time, no need to repeat that info here.
		if(bankId > 0) {
			*headerBank << addVariable(HEADER_BANK_TAG, bankId, bank.getVarType(it->first), it->second);
		}

		// storing event number in memory
		if(it->first == "evn")
		evn = it->second;

	}
	*event << headerBank;
}

// write user infos header
void evio_output :: writeUserInfoseHeader(outputContainer* output, map<string, double> data)
{
	evioDOMNodeP userHeaderBank = evioDOMNode::createEvioDOMNode(USER_HEADER_BANK_TAG, 0);
	// bank starts from 2 cause timestamp already there
	int banknum = 1;
	for(map<string, double> :: iterator it = data.begin(); it != data.end(); it++) {
		*userHeaderBank << addVariable(USER_HEADER_BANK_TAG, banknum, "d", it->second);
		banknum++;
	}

	int minNumberOfVarsToWrite = 10;
	for(unsigned i=data.size(); i<minNumberOfVarsToWrite; i++) {
		*userHeaderBank << addVariable(USER_HEADER_BANK_TAG, banknum, "d", 0.0);
		banknum++;

	}

	*event << userHeaderBank;
}


void evio_output :: writeRFSignal(outputContainer* output, FrequencySyncSignal rfsignals, gBank bank)
{
	// creating and inserting generated particles bank  >> TAG=10 NUM=0 <<
	evioDOMNodeP rfbank = evioDOMNode::createEvioDOMNode(RF_BANK_TAG, 0);

	vector<oneRFOutput> rfs = rfsignals.getOutput();
	for(unsigned i=0; i<rfs.size(); i++) {

		// each rf signal is a different container bank
		evioDOMNodeP erfsignalContainer = evioDOMNode::createEvioDOMNode(RF_BANK_TAG + i + 1, 0);

		*erfsignalContainer << addVector(RF_BANK_TAG + i + 1, bank.getVarId("id"), bank.getVarType("id"), rfs[i].getIDs());
		*erfsignalContainer << addVector(RF_BANK_TAG + i + 1, bank.getVarId("rf"), bank.getVarType("rf"), rfs[i].getValues());

		*rfbank << erfsignalContainer;

	}

	*event  << rfbank;

}

void evio_output :: writeGenerated(outputContainer* output, vector<generatedParticle> MGP, map<string, gBank> *banksMap, vector<userInforForParticle> userInfo)
{
	double MAXP             = output->gemcOpt.optMap["NGENP"].arg;
	double SAVE_ALL_MOTHERS = output->gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg ;
	int fastMCMode          = output->gemcOpt.optMap["FASTMCMODE"].arg;  // fast mc = 2 will increase prodThreshold and maxStep to 5m

	if(fastMCMode>0) SAVE_ALL_MOTHERS = 1;

	gBank bank  = getBankFromMap("generated", banksMap);
	gBank sbank = getBankFromMap("psummary", banksMap);

	vector<int> pid;
	vector<double> px;
	vector<double> py;
	vector<double> pz;
	vector<double> vx;
	vector<double> vy;
	vector<double> vz;
	vector<double> btime;
	vector<double> multiplicity;


	for(unsigned i=0; i<MAXP && i<MGP.size(); i++)
	{
		pid.push_back(MGP[i].PID);
		px.push_back(MGP[i].momentum.getX()/MeV);
		py.push_back(MGP[i].momentum.getY()/MeV);
		pz.push_back(MGP[i].momentum.getZ()/MeV);
		vx.push_back(MGP[i].vertex.getX()/mm);
		vy.push_back(MGP[i].vertex.getY()/mm);
		vz.push_back(MGP[i].vertex.getZ()/mm);
		btime.push_back(MGP[i].time);
		multiplicity.push_back(MGP[i].multiplicity);
	}

	// creating and inserting generated particles bank  >> TAG=10 NUM=0 <<
	evioDOMNodeP generatedp = evioDOMNode::createEvioDOMNode(GENERATED_PARTICLES_BANK_TAG, 0);

	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("pid"),           bank.getVarType("pid"), pid);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("px"),            bank.getVarType("px"),  px);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("py"),            bank.getVarType("py"),  py);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("pz"),            bank.getVarType("pz"),  pz);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("vx"),            bank.getVarType("vx"),  vx);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("vy"),            bank.getVarType("vy"),  vy);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("vz"),            bank.getVarType("vz"),  vz);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("time"),          bank.getVarType("time"), btime);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("multiplicity"),  bank.getVarType("multiplicity"), multiplicity);



	if(SAVE_ALL_MOTHERS)
	{
		vector<string> dname;
		vector<int>    stat;
		vector<double> etot;
		vector<double> nphe;
		vector<double> time;

		// fast MC mode, smeared and unsmeared
		vector<double> ufpx;
		vector<double> ufpy;
		vector<double> ufpz;
		vector<double> sfpx;
		vector<double> sfpy;
		vector<double> sfpz;

		int writeFastMC = 0;

		for(unsigned int i=0; i<MAXP && i<MGP.size(); i++)
		{
			if(MGP[i].fastMC.size() > 0) {
				if(MGP[i].pSum.size() != MGP[i].fastMC.size()) {
					cout << " !!! Warning: pSum and fastMC info do not match" << endl;
				} else {
					writeFastMC = 1;
				}
			}

			for(unsigned d=0; d<MGP[i].pSum.size(); d++)
			{
				dname.push_back(MGP[i].pSum[d].dname);
				stat.push_back(MGP[i].pSum[d].stat);
				etot.push_back(MGP[i].pSum[d].etot);
				nphe.push_back(MGP[i].pSum[d].nphe);
				time.push_back(MGP[i].pSum[d].t);


				if(writeFastMC) {
					ufpx.push_back(MGP[i].fastMC[d].pOrig.getX()/MeV);
					ufpy.push_back(MGP[i].fastMC[d].pOrig.getY()/MeV);
					ufpz.push_back(MGP[i].fastMC[d].pOrig.getZ()/MeV);
					sfpx.push_back(MGP[i].fastMC[d].pSmear.getX()/MeV);
					sfpy.push_back(MGP[i].fastMC[d].pSmear.getY()/MeV);
					sfpz.push_back(MGP[i].fastMC[d].pSmear.getZ()/MeV);
				}
				
			}
		}

		evioDOMNodeP summaryBank = evioDOMNode::createEvioDOMNode(GENERATED_SUMMARY_BANK_TAG, 0);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("dname"), bank.getVarType("dname"), dname);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("stat"),  bank.getVarType("stat"),  stat);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("etot"),  bank.getVarType("etot"),  etot);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("nphe"),  bank.getVarType("nphe"),  etot);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("t"),     bank.getVarType("t"),     time);

		if(writeFastMC) {
			*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("upx"),   bank.getVarType("upx"),    ufpx);
			*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("upy"),   bank.getVarType("upy"),    ufpy);
			*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("upz"),   bank.getVarType("upz"),    ufpz);
			*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("spx"),   bank.getVarType("spx"),    sfpx);
			*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("spy"),   bank.getVarType("spy"),    sfpy);
			*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("spz"),   bank.getVarType("spz"),    sfpz);
		}



		*generatedp << summaryBank;
	}

	*event << generatedp;

	// user information
	// storing 25 info / particles
	// if no user infos, these are all 0s
	evioDOMNodeP userInfoBank = evioDOMNode::createEvioDOMNode(GENERATED_USE_INFO_TAG, 0);
	for(unsigned i=0; i<25; i++) {
		// particle is the index
		vector<double> userVar(userInfo.size(), 0);

		for(unsigned p=0; p<userInfo.size(); p++) {
			if(userInfo[p].infos.size() > i) {
				userVar[p] = userInfo[p].infos[i];
			}
		}
		*userInfoBank << addVector(GENERATED_USE_INFO_TAG, i+1, "d", userVar);
	}
	*event << userInfoBank;

}

void evio_output :: initBank(outputContainer* output, gBank thisHitBank, int what)
{

	static int oldEvn;

	// new event, initialize everything
	if(oldEvn != evn)
	{
		insideBank.clear();
		insideRawIntBank.clear();
		insideDgtIntBank.clear();
		insideRawAllBank.clear();
		oldEvn = evn;
	}

	if(!insideBank[thisHitBank.bankName])
	{
		// master bank
		detectorBank = evioDOMNode::createEvioDOMNode(thisHitBank.idtag, DETECTOR_BANK_ID);
		*event << detectorBank;
		insideBank[thisHitBank.bankName] = 1;
	}

	// true information (integrated)
	if(!insideRawIntBank[thisHitBank.bankName] && what == RAWINT_ID)
	{
		detectorRawIntBank[thisHitBank.bankName] = evioDOMNode::createEvioDOMNode(thisHitBank.idtag + RAWINT_ID, 0);
		*detectorBank << detectorRawIntBank[thisHitBank.bankName];
		insideRawIntBank[thisHitBank.bankName] = 1;
	}

	// digitized information
	if(!insideDgtIntBank[thisHitBank.bankName] && what == DGTINT_ID)
	{
		detectorDgtIntBank[thisHitBank.bankName] = evioDOMNode::createEvioDOMNode(thisHitBank.idtag + DGTINT_ID, 0);
		*detectorBank << detectorDgtIntBank[thisHitBank.bankName];
		insideDgtIntBank[thisHitBank.bankName] = 1;
	}

	// true information (step by step)
	if(!insideRawAllBank[thisHitBank.bankName] && what == RAWSTEP_ID)
	{
		detectorRawAllBank[thisHitBank.bankName] = evioDOMNode::createEvioDOMNode(thisHitBank.idtag + RAWSTEP_ID, 0);
		*detectorBank << detectorRawAllBank[thisHitBank.bankName];
		insideRawAllBank[thisHitBank.bankName] = 1;
	}

	// charge / time information (step by step)
	if(!insideChargeTimeBank[thisHitBank.bankName] && what == CHARGE_TIME_ID)
	{
		detectorChargeTimeBank[thisHitBank.bankName] = evioDOMNode::createEvioDOMNode(thisHitBank.idtag + CHARGE_TIME_ID, 0);
		*detectorBank << detectorChargeTimeBank[thisHitBank.bankName];
		insideChargeTimeBank[thisHitBank.bankName] = 1;
	}


}


void evio_output :: writeG4RawIntegrated(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank rawBank = getBankFromMap("raws", banksMap);

	initBank(output, thisHitBank, RAWINT_ID);

	for(map<int, string>::iterator it =  rawBank.orderedNames.begin(); it != rawBank.orderedNames.end(); it++)
	{
		int bankId   = rawBank.getVarId(it->second);
		int bankType = rawBank.getVarBankType(it->second);

		// we only need the first hit to get the definitions
		map<string, double> raws = HO[0].getRaws();

		if(raws.find(it->second) != raws.end() && bankId > 0 && bankType == RAWINT_ID)
		{
			vector<double> thisVar;
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, double> theseRaws = HO[nh].getRaws();
				thisVar.push_back(theseRaws[it->second]);
			}
			*(detectorRawIntBank[thisHitBank.bankName]) << addVector(rawBank.idtag + thisHitBank.idtag, bankId, rawBank.getVarType(it->second), thisVar);
		}
	}
}


void evio_output :: writeG4DgtIntegrated(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;
	
	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank dgtBank = getDgtBankFromMap(hitType, banksMap);

	initBank(output, thisHitBank, DGTINT_ID);

	for(map<int, string>::iterator it =  dgtBank.orderedNames.begin(); it != dgtBank.orderedNames.end(); it++)
	{
		int bankId   = dgtBank.getVarId(it->second);
		int bankType = dgtBank.getVarBankType(it->second);

		// we only need the first hit to get the definitions
		map<string, double> dgts = HO[0].getDgtz();

		if(dgts.find(it->second) != dgts.end() && bankId > 0 && bankType == DGTINT_ID)
		{
			vector<double> thisVar;
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, double> theseDgts = HO[nh].getDgtz();
				thisVar.push_back(theseDgts[it->second]);
			}
			*(detectorDgtIntBank[thisHitBank.bankName]) << addVector(dgtBank.idtag + thisHitBank.idtag, bankId, dgtBank.getVarType(it->second), thisVar);
		}
	}
}


// index 0: hit number
// index 1: step index
// index 2: charge at electronics
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
void evio_output :: writeChargeTime(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank    = getBankFromMap(hitType, banksMap);
	gBank chargeTimeBank = getBankFromMap("chargeTime", banksMap);

	initBank(output, thisHitBank, CHARGE_TIME_ID );

	// collecting vectors from each hit into one big array
	vector<double> allHitN;
	vector<double> allStep;
	vector<double> allCharge;
	vector<double> allTime;
	vector<double> allID;

	for(unsigned int nh=0; nh<HO.size(); nh++) {

		map<int, vector<double> > thisChargeTime = HO[nh].getChargeTime();

		vector<double> thisHitN   = thisChargeTime[0];
		vector<double> thisStep   = thisChargeTime[1];
		vector<double> thisCharge = thisChargeTime[2];
		vector<double> thisTime   = thisChargeTime[3];
		vector<double> thisID     = thisChargeTime[4];


		// hit number
		if(thisHitN.size() != 1 ) {
			cout << "  !! Error: hit number should not be a vector. Bank: " << hitType << endl;
			exit(0);
		}

		for(auto h: thisHitN)   allHitN.push_back(h);
		for(auto h: thisStep)   allStep.push_back(h);
		for(auto h: thisCharge) allCharge.push_back(h);
		for(auto h: thisTime)   allTime.push_back(h);
		for(auto h: thisID)     allID.push_back(h);
	}

	// done collecting, now writing them out
	if(chargeTimeBank.getVarBankType("id") == CHARGE_TIME_ID)
	{
		// hit number
		*(detectorChargeTimeBank[thisHitBank.bankName]) << addVector(chargeTimeBank.idtag + thisHitBank.idtag, chargeTimeBank.getVarId("hitn"), chargeTimeBank.getVarType("hitn"), allHitN);
		// step index
		*(detectorChargeTimeBank[thisHitBank.bankName]) << addVector(chargeTimeBank.idtag + thisHitBank.idtag, chargeTimeBank.getVarId("stepi"), chargeTimeBank.getVarType("stepi"), allStep);
		// vector of identifiers - have to match the translation table
		*(detectorChargeTimeBank[thisHitBank.bankName]) << addVector(chargeTimeBank.idtag + thisHitBank.idtag, chargeTimeBank.getVarId("id"), chargeTimeBank.getVarType("id"), allID);
		// charge at electronics
		*(detectorChargeTimeBank[thisHitBank.bankName]) << addVector(chargeTimeBank.idtag + thisHitBank.idtag, chargeTimeBank.getVarId("q"), chargeTimeBank.getVarType("q"), allCharge);
		// time at electronics
		*(detectorChargeTimeBank[thisHitBank.bankName]) << addVector(chargeTimeBank.idtag + thisHitBank.idtag, chargeTimeBank.getVarId("t"), chargeTimeBank.getVarType("t"), allTime);
	}


}


void evio_output :: writeG4RawAll(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank allRawsBank = getBankFromMap("allraws", banksMap);

	initBank(output, thisHitBank, RAWSTEP_ID);

	for(map<int, string>::iterator it =  allRawsBank.orderedNames.begin(); it != allRawsBank.orderedNames.end(); it++)
	{
		int bankId   = allRawsBank.getVarId(it->second);
		int bankType = allRawsBank.getVarBankType(it->second);

		// we only need the first hit to get the definitions
		map<string, vector<double> > allRaws = HO[0].getAllRaws();

		if(allRaws.find(it->second) != allRaws.end() && bankId > 0 && bankType == RAWSTEP_ID)
		{
			vector<double> thisVar;
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, vector<double> > theseRaws = HO[nh].getAllRaws();

				vector<double> theseRawsSteps = theseRaws[it->second];

				for(unsigned s=0; s<theseRawsSteps.size(); s++)
				thisVar.push_back(theseRawsSteps[s]);

			}


			*(detectorRawAllBank[thisHitBank.bankName]) << addVector(allRawsBank.idtag + thisHitBank.idtag, bankId, allRawsBank.getVarType(it->second), thisVar);
		}
	}
}

void evio_output :: writeFADCMode1( map<int,  vector<hitOutput> > HO , int ev_number){

	if(HO.size() == 0) return;

	// ==== Following variables are needed for EVIO util functions PUT16, PUT31 etc
	unsigned char *b08out;
	unsigned short *b16;
	unsigned int *b32;
	unsigned long long *b64;

	// This variable will store the buffer address when crate changes
	unsigned int *buf_crate_begin; //

	// The map of hardware data, The Key is the crate/slot/channel combination,
	// and the value is a map of FADC counts as a function of sample number.
	// Note: 1st three counts represents actually the hardware identification info, i.e. crate/slot/channel, and other elements starting
	// from 3 up to nsamples+3 represent FADC counts
	map<string, map<int, int> > hardwareData;

	// map that counts how many channels are in the crate
	map<string, int> numberOfChannelsPerCrate;

	//for(unsigned int nh=0; nh<HO.size(); nh++) {
	for (std::map<int, vector<hitOutput> >::iterator it_crate=HO.begin(); it_crate!=HO.end(); ++it_crate){
		// QuantumS is a map, KEY is an FADC sample number, and value is the FADC counts
		// NOTE 1st three elements of it (KEY = 0, 1, 2) represent crate/slot/chann, and KEYs (3, 4, ... nsampes+2 ) represent FADC counts

		for( unsigned int i_hit = 0; i_hit < it_crate->second.size(); i_hit++ ){

			map<int, int> quantumS = (it_crate->second).at(i_hit).getQuantumS();

			// Let's get hardware identifiers
			string crate = fillDigits(to_string((int) quantumS[0]), "#", 5);
			string slot  = fillDigits(to_string((int) quantumS[1]), "#", 5);
			string chann = fillDigits(to_string((int) quantumS[2]), "#", 5);

			// crate-slot-channel key
			string hardwareKey = crate + "-" + slot + "-" + chann;

			// only fill hardware if it's not present already
			// the time window of a detector could be smaller than
			// the electronic time window
			// Make sure also the vector of step times is not empty,
			if(hardwareData.find(hardwareKey) == hardwareData.end() && (it_crate->second.at(i_hit).getChargeTime()[3]).size() > 0  ) {
				hardwareData[hardwareKey] = it_crate->second.at(i_hit).getQuantumS();
			} else {
				// ======== It was checked, here we have only empty hits, or hits that way off in time, e.g. hit_t = 1200ns
//                            cout<<"Warning this hardware is already filled"<<endl;
//			    cout<<"hardwareKey is "<<hardwareKey<<"      !!!"<<endl;
//                            for( int ii = 0; ii < (it_crate->second.at(i_hit).getChargeTime()[3]).size();ii++ ){
//                            
//                               cout<<"       time = " <<(it_crate->second.at(i_hit).getChargeTime()[3]).at(ii)<<"     ";
//                            }
//                            cout<<endl;
                            
		            continue;
			}

			// We should keep track of number of channels in the crate
			// With this counter, This counter will show, whether all the channels in that crate are already processed

			string CrateKey = crate;

			if( numberOfChannelsPerCrate.find(CrateKey) == numberOfChannelsPerCrate.end() ){
				numberOfChannelsPerCrate[CrateKey] = 1;
			} else {
				numberOfChannelsPerCrate[CrateKey] ++;
			}

		}
	}

	// Now we have hardware data for all crate/slot/channel combinations, and can start filling the buffer

	int oldCrate    = -1;
	int oldSlot     = -1;

	b08out = (uint8_t*) buf;

	uint32_t *nchannels;
	uint32_t *nsamples;

	evioDOMNodeP newCrate;

	int nchannelThisCrate = 0;
	int ncrates = 0;


	for(auto &hd : hardwareData) {

		vector<string> thisHardwareKey = getStringVectorFromStringWithDelimiter(hd.first, "-");

		int crate = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[0], "#", " ")));
		int slot  = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[1], "#", " ")));
		int chann = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[2], "#", " ")));

		string scrate = fillDigits(to_string(crate), "#", 5);
		string sslot  = fillDigits(to_string(slot),  "#", 5);
		string hardwareKey = scrate + "-" + sslot;
		string sCrateKey = scrate;


		//  Check, if the crate is new crate, then save the curre
		if(oldCrate != crate) {
			oldCrate = crate;
			oldSlot = -1;  // resetting slot, it's a new crate
			nchannelThisCrate  = 0;

			buf_crate_begin = (unsigned int*)b08out;
			//buf_crate_begin = (char*)b08out;
			ncrates = ncrates + 1;

			newCrate = evioDOMNode::createEvioDOMNode(crate, 0);

			// We want to check whether FADC conf data is written into evio, if not it will
			// write data and will erase corresponding crate element from the vector, for the
			// next time this data to not be written
			std::vector<int>::iterator it_crate = find(detector_crates.begin(), detector_crates.end(), crate);

			//
			if( it_crate != detector_crates.end() ){

				int confbankbanktag = 0xe10e;

				evioDOMNodeP confbank = evioDOMNode::createEvioDOMNode<string>(confbankbanktag, crate);

				int n_slotes = 19;
				int n_chann = 16;

				string conf_parms;

				conf_parms = conf_parms + "\n";
				for( int i_sl = 0; i_sl < n_slotes; i_sl++ ){

					conf_parms = conf_parms + "FADC250_SLOT " + to_string(i_sl) + "\nFADC250_NSB " + to_string(12) + "\nFADC250_NSA " + to_string(36) + "\nFADC250_ALLCH_PED ";
					for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
						conf_parms = conf_parms + " " + to_string(101.0);
					}
					conf_parms = conf_parms +"\n";
					conf_parms = conf_parms + "FADC250_ALLCH_TET ";
					for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
						conf_parms = conf_parms + " " + to_string(20);
					}
					conf_parms = conf_parms +"\n";
					conf_parms = conf_parms + "FADC250_ALLCH_GAIN ";
					for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
						conf_parms = conf_parms + " " + to_string(1.);
					}

					conf_parms = conf_parms +"\n";

				}

				*confbank<<conf_parms;
				*newCrate << confbank;
				*event << newCrate;

				detector_crates.erase(it_crate);
			}

		}

		nchannelThisCrate = nchannelThisCrate + 1;

		// every slot has its own bank
		if(oldSlot != slot) {

			oldSlot = slot;

			PUT8(slot); // slot number
			PUT32(ev_number);   // event number
			PUT64(1);   // time stamp
			nchannels = (uint32_t*) b08out; // put channels dinamically: first, save current position
			PUT32(0);   // now reserve space for channel counter
		}


		*nchannels = *nchannels + 1;

		PUT8(chann); // channel number

		nsamples = (uint32_t*) b08out; // put multi-hit dinamically: first, save current position
		PUT32(0); // now reserve space for sample counter
		for( map<int, int>::iterator it_isample = (hd.second).begin(); it_isample != (hd.second).end(); it_isample++ ) {

			// Remember 1st three elements are crate/slot/chann, therefore we want other elements hd[3], hd[4] ... hd[nsample + 3 -1]
			if( it_isample->first < 3 ){
				continue;
			}

			//cout<<"value = "<<it_isample->second<<endl;
			PUT16(abs(it_isample->second));
			*nsamples = *nsamples + 1;
		}



		// Check if all the data under this crate is processed, if yes, the
		// data should be dumped into evio
		if( nchannelThisCrate == numberOfChannelsPerCrate[sCrateKey] ){

			//int finalNumberOfWords = (b08out - (uint8_t*)buf_crate_begin + 3) / 4;
			//int finalNumberOfWords = (b08out - (uint8_t*)buf_crate_begin + 3) / 4;
			//int finalNumberOfWords = (b08out + 4 - (uint8_t*)buf_crate_begin) / 4;

			//int padding = (uint8_t*)buf_crate_begin + 4*finalNumberOfWords - b08out;

			//*newCrate << evioDOMNode::createEvioDOMNode(banktag, 0, strlen("c,i,l,N(c,Ns)"), "c,i,l,N(c,Ns)", 63, 64, buf_crate_begin, finalNumberOfWords, padding);  // Sergei wanted num to be set 0
			*newCrate << evioDOMNode::createEvioDOMNode(fadc_mode1_banktag, 0, "c,i,l,N(c,Ns)", 63, 64, buf_crate_begin, (uint32_t*)b08out);  // Sergei wanted num to be set 0
			*event << newCrate;
		}

	}

	// ======= At this moment all the FADCMode1 data is already written, so below
	// the program should iterate over remaining elements of detector_crates, and for each one
	// write FADC_conf paratmeters into evio


	if( detector_crates.size() >=1 ){

		for( vector<int>::iterator it_crate = detector_crates.begin(); it_crate != detector_crates.end(); it_crate++ ){

			int cur_crate = *it_crate;

			newCrate = evioDOMNode::createEvioDOMNode(cur_crate, 0);

			int confbanktag = 0xe10e;

			evioDOMNodeP confbank = evioDOMNode::createEvioDOMNode<string>(confbanktag, cur_crate);  // Sergei mentioned that Num should be crate number,

			int n_slotes = 19;
			int n_chann = 16;

			string conf_parms;

			conf_parms = conf_parms + "\n";
			for( int i_sl = 0; i_sl < n_slotes; i_sl++ ){

				conf_parms = conf_parms + "FADC250_SLOT " + to_string(i_sl) +  "\nFADC250_NSB " + to_string(12) + "\nFADC250_NSA " + to_string(36) + "\nFADC250_ALLCH_PED ";
				for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
					conf_parms = conf_parms + " " + to_string(101.0);
				}
				conf_parms = conf_parms +"\n";
				conf_parms = conf_parms + "FADC250_ALLCH_TET ";
				for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
					conf_parms = conf_parms + " " + to_string(20);
				}
				conf_parms = conf_parms +"\n";
				conf_parms = conf_parms + "FADC250_ALLCH_GAIN ";
				for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
					conf_parms = conf_parms + " " + to_string(1.);
				}

				conf_parms = conf_parms +"\n";

			}

			*confbank<<conf_parms;
			*newCrate << confbank;
			*event << newCrate;
		}
		detector_crates.clear();
	}


}


// write fadc mode 1 (full signal shape) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
// This function takes as an argument vector of hitOutputs, and writes all hits into evio in a Mode1 format
void evio_output :: writeFADCMode1(outputContainer* output, vector<hitOutput> HO, int ev_number)
{
	if(HO.size() == 0) return;

	// ==== Following variables are needed for EVIO util functions PUT16, PUT31 etc
	unsigned char *b08out;
	unsigned short *b16;
	unsigned int *b32;
	unsigned long long *b64;


	// This variable will store the buffer address when crate changes
	unsigned int *buf_crate_begin; //
	//char *buf_crate_begin; //

	// The FADC Mode1 Bank tag
	int banktag = 0xe101;


	// The map of hardware data, The Key is the crate/slot/channel combination,
	// and the value is a map of FADC counts as a function of sample number.
	// Note: 1st three counts represents actually the hardware identification info, i.e. crate/slot/channel, and other elements starting
	// from 3 up to nsamples+3 represent FADC counts
	map<string, map<int, int> > hardwareData;

	// map that counts how many channels are in the crate
	map<string, int> numberOfChannelsPerCrate;

	for(unsigned int nh=0; nh<HO.size(); nh++) {

		// QuantumS is a map, KEY is an FADC sample number, and value is the FADC counts
		// NOTE 1st three elements of it (KEY = 0, 1, 2) represent crate/slot/chann, and KEYs (3, 4, ... nsampes+2 ) represent FADC counts
		map<int, int> quantumS = HO[nh].getQuantumS();

		// Let's get hardware identifiers
		string crate = fillDigits(to_string((int) quantumS[0]), "#", 5);
		string slot  = fillDigits(to_string((int) quantumS[1]), "#", 5);
		string chann = fillDigits(to_string((int) quantumS[2]), "#", 5);

		// crate-slot-channel key
		string hardwareKey = crate + "-" + slot + "-" + chann;

		// only fill hardware if it's not present already
		// the time window of a detector could be smaller than
		// the electronic time window
		// Make sure also the vector of step times is not empty,
		if(hardwareData.find(hardwareKey) == hardwareData.end() && (HO[nh].getChargeTime()[3]).size() > 0  ) {
			hardwareData[hardwareKey] = HO[nh].getQuantumS();

			vector<double> stepTimes   = HO[nh].getChargeTime()[3];

		} else {

			// ======== It was checked, here we have only empty hits, or hits that way off in time, e.g. hit_t = 1200ns
			//cout<<"hardwareKey is "<<hardwareKey<<"      !!!"<<endl;
			continue;
		}

		// We should keep track of number of channels in the crate
		// With this counter, This counter will show, whether all the channels in that crate are already processed

		string CrateKey = crate;

		if( numberOfChannelsPerCrate.find(CrateKey) == numberOfChannelsPerCrate.end() ){
			numberOfChannelsPerCrate[CrateKey] = 1;
		} else {
			numberOfChannelsPerCrate[CrateKey] ++;
		}


	}


	// Now we have hardware data for all crate/slot/channel combinations, and can start filling the buffer

	int oldCrate    = -1;
	int oldSlot     = -1;

	b08out = (uint8_t*) buf;

	uint32_t *nchannels;
	uint32_t *nsamples;

	evioDOMNodeP newCrate;

	int nchannelThisCrate = 0;
	int ncrates = 0;

	for(auto &hd : hardwareData) {

		vector<string> thisHardwareKey = getStringVectorFromStringWithDelimiter(hd.first, "-");

		int crate = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[0], "#", " ")));
		int slot  = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[1], "#", " ")));
		int chann = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[2], "#", " ")));

		string scrate = fillDigits(to_string(crate), "#", 5);
		string sslot  = fillDigits(to_string(slot),  "#", 5);
		string hardwareKey = scrate + "-" + sslot;
		string sCrateKey = scrate;


		//  Check, if the crate is new crate, then save the curre
		if(oldCrate != crate) {
			oldCrate = crate;
			oldSlot = -1;  // resetting slot, it's a new crate
			nchannelThisCrate  = 0;

			buf_crate_begin = (unsigned int*)b08out;
			//buf_crate_begin = (char*)b08out;
			ncrates = ncrates + 1;

			newCrate = evioDOMNode::createEvioDOMNode(crate, 0);

			// We want to check whether FADC conf data is written into evio, if not it will
			// write data and will erase corresponding crate element from the vector, for the
			// next time this data to not be written
			std::vector<int>::iterator it_crate = find(detector_crates.begin(), detector_crates.end(), crate);

			//
			if( it_crate != detector_crates.end() ){

				int confbankbanktag = 0xe10e;

				evioDOMNodeP confbank = evioDOMNode::createEvioDOMNode<string>(confbankbanktag, crate);

				int n_slotes = 19;
				int n_chann = 16;

				string conf_parms;

				conf_parms = conf_parms + "\n";
				for( int i_sl = 0; i_sl < n_slotes; i_sl++ ){

					conf_parms = conf_parms + "FADC250_SLOT " + to_string(i_sl) + "\nFADC250_NSB " + to_string(12) + "\nFADC250_NSA " + to_string(36) + "\nFADC250_ALLCH_PED ";
					for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
						conf_parms = conf_parms + " " + to_string(101.0);
					}
					conf_parms = conf_parms +"\n";
					conf_parms = conf_parms + "FADC250_ALLCH_TET ";
					for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
						conf_parms = conf_parms + " " + to_string(20);
					}
					conf_parms = conf_parms +"\n";
					conf_parms = conf_parms + "FADC250_ALLCH_GAIN ";
					for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
						conf_parms = conf_parms + " " + to_string(1.);
					}

					conf_parms = conf_parms +"\n";

				}

				*confbank<<conf_parms;
				*newCrate << confbank;
				*event << newCrate;

				detector_crates.erase(it_crate);
			}

		}

		nchannelThisCrate = nchannelThisCrate + 1;

		// every slot has its own bank
		if(oldSlot != slot) {

			oldSlot = slot;

			PUT8(slot); // slot number
			PUT32(ev_number);   // event number
			PUT64(1);   // time stamp
			nchannels = (uint32_t*) b08out; // put channels dinamically: first, save current position
			PUT32(0);   // now reserve space for channel counter
		}


		*nchannels = *nchannels + 1;

		PUT8(chann); // channel number

		nsamples = (uint32_t*) b08out; // put multi-hit dinamically: first, save current position
		PUT32(0); // now reserve space for sample counter
		for( map<int, int>::iterator it_isample = (hd.second).begin(); it_isample != (hd.second).end(); it_isample++ ) {

			// Remember 1st three elements are crate/slot/chann, therefore we want other elements hd[3], hd[4] ... hd[nsample + 3 -1]
			if( it_isample->first < 3 ){
				continue;
			}

			//cout<<"value = "<<it_isample->second<<endl;
			PUT16(abs(it_isample->second));
			*nsamples = *nsamples + 1;
		}



		// Check if all the data under this crate is processed, if yes, the
		// data should be dumped into evio
		if( nchannelThisCrate == numberOfChannelsPerCrate[sCrateKey] ){

			//int finalNumberOfWords = (b08out - (uint8_t*)buf_crate_begin + 3) / 4;
			//int finalNumberOfWords = (b08out - (uint8_t*)buf_crate_begin + 3) / 4;
			//int finalNumberOfWords = (b08out + 4 - (uint8_t*)buf_crate_begin) / 4;

			//int padding = (uint8_t*)buf_crate_begin + 4*finalNumberOfWords - b08out;

			//*newCrate << evioDOMNode::createEvioDOMNode(banktag, 0, strlen("c,i,l,N(c,Ns)"), "c,i,l,N(c,Ns)", 63, 64, buf_crate_begin, finalNumberOfWords, padding);  // Sergei wanted num to be set 0
			*newCrate << evioDOMNode::createEvioDOMNode(banktag, 0, "c,i,l,N(c,Ns)", 63, 64, buf_crate_begin, (uint32_t*)b08out);  // Sergei wanted num to be set 0
			*event << newCrate;
		}

	}

	// ======= At this moment all the FADCMode1 data is already written, so below
	// the program should iterate over remaining elements of detector_crates, and for each one
	// write FADC_conf paratmeters into evio


	if( detector_crates.size() >=1 ){

		for( vector<int>::iterator it_crate = detector_crates.begin(); it_crate != detector_crates.end(); it_crate++ ){

			int cur_crate = *it_crate;

			newCrate = evioDOMNode::createEvioDOMNode(cur_crate, 0);

			int confbanktag = 0xe10e;

			evioDOMNodeP confbank = evioDOMNode::createEvioDOMNode<string>(confbanktag, cur_crate);  // Sergei mentioned that Num should be crate number,

			int n_slotes = 19;
			int n_chann = 16;

			string conf_parms;

			conf_parms = conf_parms + "\n";
			for( int i_sl = 0; i_sl < n_slotes; i_sl++ ){

				conf_parms = conf_parms + "FADC250_SLOT " + to_string(i_sl) +  "\nFADC250_NSB " + to_string(12) + "\nFADC250_NSA " + to_string(36) + "\nFADC250_ALLCH_PED ";
				for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
					conf_parms = conf_parms + " " + to_string(101.0);
				}
				conf_parms = conf_parms +"\n";
				conf_parms = conf_parms + "FADC250_ALLCH_TET ";
				for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
					conf_parms = conf_parms + " " + to_string(20);
				}
				conf_parms = conf_parms +"\n";
				conf_parms = conf_parms + "FADC250_ALLCH_GAIN ";
				for( int i_ch = 0; i_ch < n_chann; i_ch++ ){
					conf_parms = conf_parms + " " + to_string(1.);
				}

				conf_parms = conf_parms +"\n";

			}

			*confbank<<conf_parms;
			*newCrate << confbank;
			*event << newCrate;
		}
		detector_crates.clear();
	}


}


// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void evio_output :: writeFADCMode7(outputContainer* output, vector<hitOutput> HO, int ev_number)
{
	unsigned char *b08out;
	unsigned short *b16;
	unsigned int *b32;
	unsigned long long *b64;

	int banktag = 0xe102;

	map<string, map<int, int> > hardwareData;

	// map that counts how many channels are active / slot
	map<string, int> numberOfChannelsPerSlot;

	// first, reorder all hits into a map
	for(unsigned int nh=0; nh<HO.size(); nh++) {

		vector<double> hardware = HO[nh].getChargeTime()[5];
		vector<double> identifier = HO[nh].getChargeTime()[4];

		string crate = fillDigits(to_string((int) hardware[0]), "#", 5);
		string slot  = fillDigits(to_string((int) hardware[1]), "#", 5);
		string chann = fillDigits(to_string((int) hardware[2]), "#", 5);

		// crate-slot-channel key
		string hardwareKey = crate + "-" + slot + "-" + chann;

		// only fill hardware if it's not present already
		// the time window of a detector could be smaller than
		// the electronic time window
		if(hardwareData.find(hardwareKey) == hardwareData.end()) {
			hardwareData[hardwareKey] = HO[nh].getQuantumS();
		} else {
			continue;
		}

		// crate-slot key
		hardwareKey = crate + "-" + slot;

		// keeping track of how many channels / slot so we know when to write data
		if(numberOfChannelsPerSlot.find(hardwareKey) == numberOfChannelsPerSlot.end()) {
			numberOfChannelsPerSlot[hardwareKey] = 1;
		} else {
			numberOfChannelsPerSlot[hardwareKey]++;
		}
	}

	int oldCrate    = -1;
	int oldSlot     = -1;

	b08out = (uint8_t*) buf;

	uint32_t *nchannels;
	uint32_t *nhits;

	evioDOMNodeP newCrate;

	int nchannelThisSlot = 0;

	for(auto &hd : hardwareData) {

		vector<string> thisHardwareKey = getStringVectorFromStringWithDelimiter(hd.first, "-");

		int crate = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[0], "#", " ")));
		int slot  = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[1], "#", " ")));
		int chann = stoi(trimSpacesFromString(replaceCharInStringWithChars(thisHardwareKey[2], "#", " ")));

		string scrate = fillDigits(to_string(crate), "#", 5);
		string sslot  = fillDigits(to_string(slot),  "#", 5);
		string hardwareKey = scrate + "-" + sslot;

		if(oldCrate != crate) {
			oldCrate = crate;
			oldSlot = -1;  // resetting slot, it's a new crate

			//cout << " creating new crate node " << endl;

			newCrate = evioDOMNode::createEvioDOMNode(crate, 1);
		}

		// every slot has its own bank
		if(oldSlot != slot) {

			nchannelThisSlot = 0;

			// cout << " !  crate " << crate << "  slot " << slot << "  channel " << chann << " totChannels " << numberOfChannelsPerSlot[hardwareKey] << " count so far " << nchannelThisSlot << " " << hardwareKey << endl;

			oldSlot = slot;
			
			PUT8(slot); // slot number
			PUT32(ev_number);   // event number
			PUT64(1);   // time stamp
			nchannels = (uint32_t*) b08out; // put channels dinamically: first, save current position
			PUT32(0);   // now reserve space for channel counter

		}

		// every entry is a new channel
		nchannelThisSlot++;

		*nchannels = *nchannels + 1;

		PUT8(chann); // channel number
		nhits = (uint32_t*) b08out; // put multi-hit dinamically: first, save current position
		PUT32(0); // now reserve space for hit counter
		*nhits = *nhits + 1; // in our case, only 1 hit
		PUT16(2); // pulse time
		PUT32(3); // pulse integral
		PUT16(4); // pulse min
		PUT16(5); // pulse max


		// cout << " >  crate " << crate << "  slot " << slot << "  channel " << chann << " totChannels " << numberOfChannelsPerSlot[hardwareKey] << " count so far " << nchannelThisSlot << " " << hardwareKey  << endl;

		// channel is new, writing
		if(nchannelThisSlot == numberOfChannelsPerSlot[hardwareKey]) {

			int finalNumberOfWords = (b08out - (uint8_t*)buf + 3) / 4;

			// filling crate bank
			*newCrate << evioDOMNode::createEvioDOMNode(banktag, 0, 65, "c,i,l,N(c,N(s,i,s,s))", 66, 67, buf,finalNumberOfWords); // Sergei wanted num to be set 0
			*event << newCrate;

		}
	}
}




evioDOMNodeP addVariable(int tag, int num, string type, double value)
{
	// return right away if "d"
	if(type == "d")
	return evioDOMNode::createEvioDOMNode(tag, num, &value, 1);

	// otherwise check
	if(type == "i")
	{
		int    varI = (int) value;
		return evioDOMNode::createEvioDOMNode(tag, num, &varI, 1);
	}

	// repeated return to make compiler happy
	return evioDOMNode::createEvioDOMNode(tag, num, &value, 1);

}


evioDOMNodeP addVariable(int tag, int num, string type, int value)
{
	return evioDOMNode::createEvioDOMNode(tag, num, &value, 1);
}

evioDOMNodeP addVariable(int tag, int num, string type, string value)
{
	return evioDOMNode::createEvioDOMNode(tag, num, &value, 1);
}


evioDOMNodeP addVector(int tag, int num, string type, vector<double> value)
{
	// return right away if "d"
	if(type == "d")
	return evioDOMNode::createEvioDOMNode(tag, num, value);

	// otherwise convert the double into int
	if(type == "i")
	{
		vector<int> VI;
		for(unsigned i=0; i<value.size(); i++)
		VI.push_back(value[i]);

		return evioDOMNode::createEvioDOMNode(tag, num, VI);
	}

	return evioDOMNode::createEvioDOMNode(tag, num, value);
}

// we don't convert back to double, that should be in the
// variable definition
evioDOMNodeP addVector(int tag, int num, string type, vector<int> value)
{
	return evioDOMNode::createEvioDOMNode(tag, num, value);
}


evioDOMNodeP addVector(int tag, int num, string type, vector<string> value)
{
	return evioDOMNode::createEvioDOMNode(tag, num, value);
}




void evio_output :: writeEvent(outputContainer* output)
{
	output->pchan->write(*event);
	delete event;
}

