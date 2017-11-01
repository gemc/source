// gemc headers
#include "txt_simple_output.h"
#include "utils.h"

// C++ headers
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

#define INDENT_DEPTH 4


string indent(int depth)
{
	return string(depth * INDENT_DEPTH, ' ');
}


// record the simulation conditions
// the format is a string for each variable
void txt_simple_output :: recordSimConditions(outputContainer* output, map<string, string> simcons)
{
	ofstream *txtout = output->txtoutput ;

	*txtout << "Simulation Conditions, TAG " << SIMULATION_CONDITIONS_BANK_TAG << " {" << endl;

	for(map<string, string>::iterator it = simcons.begin(); it != simcons.end(); it++)
		*txtout << indent(1) << it->first << " " << it->second << endl;

	*txtout << "}" << endl;
}




// write header bank
// initialize insideBank to zero
void txt_simple_output :: writeHeader(outputContainer* output, map<string, double> data, gBank bank)
{
	insideBank.clear();
	ofstream *txtout = output->txtoutput ;

	*txtout << "Event {" << endl;
	*txtout << indent(1) << "Header Bank {" << endl;
	*txtout << indent(2) << "(" << HEADER_BANK_TAG << ", " << bank.getVarId("time") << ") " << "time:\t" << timeStamp() << endl;

	for(map<string, double> :: iterator it = data.begin(); it != data.end(); it++)
	{
		int bankId = bank.getVarId(it->first);

		if(bankId > 0)
			*txtout << indent(2) << "(" << HEADER_BANK_TAG << ", " << bankId << ") " << it->first << ":\t" << it->second << endl;

	}
	*txtout << indent(1) << "}" << endl;
}


// write user header bank
// initialize insideBank to zero
void txt_simple_output :: writeUserInfoseHeader(outputContainer* output, map<string, double> data)
{

	ofstream *txtout = output->txtoutput ;

	*txtout << indent(1) << "User Header Bank {" << endl;

	int banknum = 1;
	for(map<string, double> :: iterator it = data.begin(); it != data.end(); it++) {
		*txtout << indent(2) << "(" << USER_HEADER_BANK_TAG << ", " << banknum << ") " << it->first << ":\t" << it->second << endl;
		banknum++;
	}
	*txtout << indent(1) << "}" << endl;
}




void txt_simple_output :: writeRFSignal(outputContainer* output, FrequencySyncSignal rfsignals, gBank bank)
{
	ofstream *txtout = output->txtoutput ;

	*txtout << indent(1) << "RF Signals Bank {" << endl;


	vector<oneRFOutput> rfs = rfsignals.getOutput();

	for(unsigned int i=0; i<rfs.size();  i++) {

		vector<int>    rfid = rfs[i].getIDs();
		vector<double> rfva = rfs[i].getValues();

		for(unsigned int j=0; j<rfid.size();  j++) {
			*txtout << indent(2) << "RF " << i+1 << " ID: "    << rfid[j] << "   Signal: " << rfva[j] << endl;
		}


	}
	*txtout << indent(1) << "}" << endl;
}

string plural_or_sing(string stem, int n){
	if (n == 1)
		return stem + "";
	return stem + "s";
}


void txt_simple_output :: writeGenerated(outputContainer* output, vector<generatedParticle> MGP, map<string, gBank> *banksMap, vector<userInforForParticle> userInfo)
{
	double MAXP = output->gemcOpt.optMap["NGENP"].arg;
	ofstream *txtout = output->txtoutput ;

	gBank bank = getBankFromMap("generated", banksMap);

	*txtout << indent(1) << "Generated Particles Bank {" << endl;

	int writeFastMC = 0;

	for(unsigned int i=0; i<MAXP && i<MGP.size(); i++)
	{
		*txtout << indent(2) << "particle " << i+1 <<" {" << endl;
		*txtout << indent(3) << "pid:\t" << MGP[i].PID << endl;
		*txtout << indent(3) << "mom:\t" << MGP[i].momentum/MeV << " MeV" << endl;
		*txtout << indent(3) << "vert:\t" << MGP[i].vertex/mm << " mm" << endl;

		if(MGP[i].fastMC.size() > 0) {
			if(MGP[i].pSum.size() != MGP[i].fastMC.size()) {
				cout << " !!! Warning: pSum and fastMC info do not match" << endl;
			} else {
				writeFastMC = 1;
			}
		}

		for(unsigned d=0; d<MGP[i].pSum.size(); d++)
		{
			*txtout << indent(3) << "generated {" << endl;
			*txtout << indent(4) << "hits:\t" << MGP[i].pSum[d].stat << endl;
			*txtout << indent(4) << "dname:\t" << MGP[i].pSum[d].dname << endl;

			if(MGP[i].pSum[d].etot > 0)
			{
				*txtout << indent(4) << "etot:\t" << MGP[i].pSum[d].etot/MeV << " MeV" << endl;
			    *txtout << indent(4) << "time:\t" << MGP[i].pSum[d].t << " ns" << endl;
			}
			else if(MGP[i].pSum[d].nphe > 0)
			{
				*txtout << indent(4) << "nphe:\t" << MGP[i].pSum[d].nphe << " nphe" << endl;
			    *txtout << indent(4) << "time:\t" << MGP[i].pSum[d].t << " ns" << endl;
			}
			if(writeFastMC == 1) {
				*txtout << indent(4) << "orig mom:\t" << MGP[i].fastMC[d].pOrig/MeV  << " MeV" << endl;
				*txtout << indent(4) << "smeared mom:\t(" << MGP[i].fastMC[d].pSmear.x() << ", "
													     <<  MGP[i].fastMC[d].pSmear.y() << ", "
														 << MGP[i].fastMC[d].pSmear.z() <<
														 ") GeV" << endl;

			}
			*txtout << indent(3) << "}" << endl;
		}
		*txtout << indent(2) << "}" << endl;
	}

	if(userInfo.size()) {
		*txtout << indent(2) << "Generated User Info Particles Bank:" << endl;
		for(unsigned p=0; p<userInfo.size(); p++) {
			*txtout << "      ";
			for(unsigned i=0; i<userInfo[p].infos.size(); i++) {
				*txtout << userInfo[p].infos[i] << " " ;
			}
			*txtout << indent(2) << endl;
		}
	}


	*txtout << indent(1) << "}" << endl;
}

void txt_simple_output ::  initBank(outputContainer* output, gBank thisHitBank)
{
}

// write out true information. This is common to all banks
// and not contained in the banks definitions
void txt_simple_output ::  writeG4RawIntegrated(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank rawBank = getBankFromMap("raws", banksMap);

	ofstream *txtout = output->txtoutput ;
	*txtout << indent(1) << thisHitBank.bankName << " (" << thisHitBank.idtag << ", " << DETECTOR_BANK_ID << ") integrated true infos bank (" << thisHitBank.idtag + RAWINT_ID << ", 0) {" << endl;

	for(map<int, string>::iterator it =  rawBank.orderedNames.begin(); it != rawBank.orderedNames.end(); it++)
	{
		int bankId   = rawBank.getVarId(it->second);
		int bankType = rawBank.getVarBankType(it->second);

		// we only need the first hit to get the definitions
		map<string, double> raws = HO[0].getRaws();

		// bankID 0 is hit index
		if(raws.find(it->second) != raws.end() && bankId >= 0 && bankType == RAWINT_ID)
		{
			*txtout << indent(2) << "(" << rawBank.idtag + thisHitBank.idtag << ", " << bankId << ") " << it->second << ":\t" ;
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, double> theseRaws = HO[nh].getRaws();
				*txtout <<  std::setprecision(12) << theseRaws[it->second] << "\t" ;
			}
			*txtout << endl;
		}
	}

	*txtout << indent(1) << "}" << endl;
}




// write out true information step by step. This is common to all banks
// and not contained in the banks definitions
void txt_simple_output ::  writeG4RawAll(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank allRawsBank = getBankFromMap("allraws", banksMap);

	ofstream *txtout = output->txtoutput ;
	*txtout << indent(1) << thisHitBank.bankName << " (" << thisHitBank.idtag << ", " << DETECTOR_BANK_ID << ") step by step true infos bank (" << thisHitBank.idtag + RAWSTEP_ID << ", 0) {" << endl;

	for(map<int, string>::iterator it =  allRawsBank.orderedNames.begin(); it != allRawsBank.orderedNames.end(); it++)
	{
		int bankId   = allRawsBank.getVarId(it->second);
		int bankType = allRawsBank.getVarBankType(it->second);

		// we only need the first hit to get the definitions
		map<string, vector<double> > allRaws = HO[0].getAllRaws();

		// bankID 0 is hit index
		if(allRaws.find(it->second) != allRaws.end() && bankId >= 0 && bankType == RAWSTEP_ID)
		{
			*txtout << indent(2) << "(" << allRawsBank.idtag + thisHitBank.idtag << ", " << bankId << ") " << it->second << ":\t" ;
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, vector<double> > theseRaws = HO[nh].getAllRaws();

				vector<double> theseRawsSteps = theseRaws[it->second];

				for(unsigned s=0; s<theseRawsSteps.size(); s++)
					*txtout << theseRawsSteps[s] << "\t" ;
			}
			*txtout << endl;
		}
	}
	*txtout << indent(1) << "}" << endl;
}


void txt_simple_output ::  writeG4DgtIntegrated(outputContainer* output, vector<hitOutput> HO,  string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank dgtBank = getDgtBankFromMap(hitType, banksMap);

	ofstream *txtout = output->txtoutput ;
	*txtout << indent(1) << thisHitBank.bankName << " (" << thisHitBank.idtag << ", " << DETECTOR_BANK_ID << ") integrated digitized bank (" << thisHitBank.idtag + DGTINT_ID << ", 0) {" << endl;

	for(map<int, string>::iterator it =  dgtBank.orderedNames.begin(); it != dgtBank.orderedNames.end(); it++)
	{
		int bankId   = dgtBank.getVarId(it->second);
		int bankType = dgtBank.getVarBankType(it->second);

		// we only need the first hit to get the definitions
		map<string, double> dgts = HO[0].getDgtz();

		// bankID 0 is hit index
		if(dgts.find(it->second) != dgts.end() && bankId > 0 && bankType == DGTINT_ID)
		{
			*txtout << indent(2) << "(" << dgtBank.idtag + thisHitBank.idtag << ", " << bankId << ") " << it->second << ":\t";

			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, double> theseDgts = HO[nh].getDgtz();
				*txtout << theseDgts[it->second] << "\t" ;
			}
			*txtout << endl;
		}
	}

	*txtout << indent(1) << "}" << endl;
}

// index 0: hit number
// index 1: step index
// index 2: charge at electronics
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
void txt_simple_output :: writeChargeTime(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	if(HO.size() == 0) return;

	gBank thisHitBank    = getBankFromMap(hitType, banksMap);
	gBank chargeTimeBank = getBankFromMap("chargeTime", banksMap);

	ofstream *txtout = output->txtoutput ;
	*txtout << indent(1) << thisHitBank.bankName << " (" << thisHitBank.idtag << ", " << DETECTOR_BANK_ID << ") charge time infos (as seen by electronics) bank  (" << thisHitBank.idtag + CHARGE_TIME_ID << ", 0) {" << endl;

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
		*txtout << indent (2) << "hit " << thisHitN[0] << " {" << endl;
		*txtout << indent (3) << "steps:\t" << thisStep.size() << endl;

		// identifier
		*txtout << indent (3) << "Identifier:\t " ;
		for(unsigned s=0; s<thisID.size(); s++) *txtout << (int) thisID[s] << " "  ;
		*txtout << endl ;


		// step index
		*txtout << indent (3) << "Step number:\t " ;
		for(unsigned s=0; s<thisStep.size(); s++) *txtout << (int) thisStep[s] << " "  ;
		*txtout << endl ;

		// charge at electronics
		*txtout << indent (3) << "Charge at electronics:\t " ;
		for(unsigned s=0; s<thisCharge.size(); s++) *txtout << thisCharge[s] << " "  ;
		*txtout << endl ;

		// time at electronics
		*txtout << indent (3) << "Time at electronics:\t " ;
		for(unsigned s=0; s<thisTime.size(); s++) *txtout << thisTime[s] << " "  ;
		*txtout << endl ;

		*txtout << indent (2) << "}" << endl;

	}
	*txtout << indent(1) << "}" << endl;

}

// write fadc mode 1 (full signal shape) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void txt_simple_output :: writeFADCMode1(outputContainer* output, vector<hitOutput> HO, int event_number)
{
}


void txt_simple_output :: writeFADCMode1( map<int, vector<hitOutput> >, int)
{
}

// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void txt_simple_output :: writeFADCMode7(outputContainer* output, vector<hitOutput> HO, int event_number)
{
}


void txt_simple_output :: writeEvent(outputContainer* output)
{
	ofstream *txtout = output->txtoutput ;
	*txtout << "}" << endl;
}
