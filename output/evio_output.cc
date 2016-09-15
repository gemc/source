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

	for(map<string, double> :: iterator it = data.begin(); it != data.end(); it++)
	{
		int bankId = bank.getVarId(it->first);

		if(bankId)
		*headerBank << addVariable(HEADER_BANK_TAG, bankId, bank.getVarType(it->first), it->second);

		// storing event number in memory
		if(it->first == "evn")
		evn = it->second;

	}
	*event << headerBank;
}

void evio_output :: writeRFSignal(outputContainer* output, FrequencySyncSignal rfsignals, gBank bank)
{
	// creating and inserting generated particles bank  >> TAG=10 NUM=0 <<
	evioDOMNodeP rfbank = evioDOMNode::createEvioDOMNode(RF_BANK_TAG, 0);

	vector<oneRFOutput> rfs = rfsignals.getOutput();
	for(unsigned i=0; i<rfs.size(); i++) {

		// each rf signal is a different bank
		evioDOMNodeP erfsignal = evioDOMNode::createEvioDOMNode(RF_BANK_TAG, i+1);

		*erfsignal << addVector(RF_BANK_TAG, bank.getVarId("id"), bank.getVarType("id"), rfs[i].getIDs());
		*erfsignal << addVector(RF_BANK_TAG, bank.getVarId("rf"), bank.getVarType("rf"), rfs[i].getValues());

		*rfbank << erfsignal;

	}

	*event  << rfbank;

}

void evio_output :: writeGenerated(outputContainer* output, vector<generatedParticle> MGP, map<string, gBank> *banksMap)
{
	double MAXP             = output->gemcOpt.optMap["NGENP"].arg;
	double SAVE_ALL_MOTHERS = output->gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg ;

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

		for(unsigned int i=0; i<MAXP && i<MGP.size(); i++)
		{


			for(unsigned d=0; d<MGP[i].pSum.size(); d++)
			{
				dname.push_back(MGP[i].pSum[d].dname);
				stat.push_back(MGP[i].pSum[d].stat);
				etot.push_back(MGP[i].pSum[d].etot);
				nphe.push_back(MGP[i].pSum[d].nphe);
				time.push_back(MGP[i].pSum[d].t);
			}

		}

		evioDOMNodeP summaryBank = evioDOMNode::createEvioDOMNode(GENERATED_SUMMARY_BANK_TAG, 0);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("dname"), bank.getVarType("dname"), dname);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("stat"),  bank.getVarType("stat"),  stat);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("etot"),  bank.getVarType("etot"),  etot);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("nphe"),  bank.getVarType("nphe"),  etot);
		*summaryBank << addVector(GENERATED_SUMMARY_BANK_TAG, sbank.getVarId("t"),     bank.getVarType("t"),     time);
		*generatedp << summaryBank;
	}

	*event << generatedp;

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

// write fadc mode 1 (full signal shape) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void evio_output :: writeFADCMode1(outputContainer* output, vector<hitOutput> HO)
{
}


// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void evio_output :: writeFADCMode7(outputContainer* output, vector<hitOutput> HO)
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
			PUT32(1);   // event number
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
			*newCrate << evioDOMNode::createEvioDOMNode(banktag, crate, 65, "c,i,l,N(c,N(s,i,s,s))", 66, 67, buf,finalNumberOfWords);
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

