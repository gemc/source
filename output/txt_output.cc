// gemc headers
#include "txt_output.h"
#include "utils.h"

// C++ headers
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


// record the simulation conditions
// the format is a string for each variable
void txt_output :: recordSimConditions(outputContainer* output, map<string, string> simcons)
{
	ofstream *txtout = output->txtoutput ;
	
	*txtout << "   Simulation Conditions, TAG " << SIMULATION_CONDITIONS_BANK_TAG << ":" << endl;
	
	for(map<string, string>::iterator it = simcons.begin(); it != simcons.end(); it++)
		*txtout << "   > " << it->first << " " << it->second << endl;
}




// write header bank
// initialize insideBank to zero
void txt_output :: writeHeader(outputContainer* output, map<string, double> data, gBank bank)
{
	insideBank.clear();
	ofstream *txtout = output->txtoutput ;
	
	
	*txtout << " --- Header Bank -- " << endl;
	*txtout << "    - (" << HEADER_BANK_TAG << ", " << bank.getVarId("time") << ") " << "time:\t" << timeStamp() << endl;
	
	for(map<string, double> :: iterator it = data.begin(); it != data.end(); it++)
	{
		int bankId = bank.getVarId(it->first);
		
		if(bankId)
			*txtout << "    - (" << HEADER_BANK_TAG << ", " << bankId << ") " << it->first << ":\t" << it->second << endl;
	}
	*txtout << " --- End of Header Bank -- " << endl;
	
	
}

void txt_output :: writeRFSignal(outputContainer* output, FrequencySyncSignal rfsignals, gBank bank)
{
	ofstream *txtout = output->txtoutput ;

	*txtout << " --- RF Signals Bank -- " << endl;


	vector<oneRFOutput> rfs = rfsignals.getOutput();

	for(unsigned int i=0; i<rfs.size();  i++) {

		vector<int>    rfid = rfs[i].getIDs();
		vector<double> rfva = rfs[i].getValues();

		for(unsigned int j=0; j<rfid.size();  j++) {
			*txtout << "    - RF " << i+1 << " ID: "    << rfid[j] << "   Signal: " << rfva[j] << endl;
		}


	}
	*txtout << " --- End of RF Signals Bank Bank -- " << endl;
}

void txt_output :: writeGenerated(outputContainer* output, vector<generatedParticle> MGP, map<string, gBank> *banksMap)
{
	double MAXP = output->gemcOpt.optMap["NGENP"].arg;
	ofstream *txtout = output->txtoutput ;

	gBank bank = getBankFromMap("generated", banksMap);

	*txtout << " --- Generated Particles Bank -- " << endl;
	
	for(unsigned int i=0; i<MAXP && i<MGP.size(); i++)
	{
		*txtout << "    - Particle "    << i+1 << " pid: "<< MGP[i].PID
		        << "   -  mom: "        << MGP[i].momentum/MeV
		        << " MeV   -  vert: "   << MGP[i].vertex/mm << " mm" << endl;
				
		for(unsigned d=0; d<MGP[i].pSum.size(); d++)
		{
			*txtout << "    - Hit >"     << MGP[i].pSum[d].dname
			        << "< Has "          << MGP[i].pSum[d].stat << " hit";
			if(MGP[i].pSum[d].stat>1) *txtout << "s";
			
			if(MGP[i].pSum[d].etot > 0)
			{
				*txtout << " with etot "     << MGP[i].pSum[d].etot/MeV
			            << " MeV and time "  << MGP[i].pSum[d].t << " ns" << endl;
			}
			else if(MGP[i].pSum[d].nphe > 0)
			{
				*txtout << " with nphe "      << MGP[i].pSum[d].nphe
			            << " nphe and time "  << MGP[i].pSum[d].t << " ns" << endl;
			}
		}
	}
	*txtout << " --- End of Generated Particles Bank -- " << endl;
}

void txt_output ::  initBank(outputContainer* output, gBank thisHitBank)
{
	if(!insideBank[thisHitBank.bankName])
	{
		ofstream *txtout = output->txtoutput ;
		*txtout << " --- " << thisHitBank.bankName << "  (" << thisHitBank.idtag << ", " << DETECTOR_BANK_ID << ") ---- " << endl;
		insideBank[thisHitBank.bankName] = 1;
	}
}

// write out true information. This is common to all banks
// and not contained in the banks definitions
void txt_output ::  writeG4RawIntegrated(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank rawBank = getBankFromMap("raws", banksMap);

	initBank(output, thisHitBank);
	ofstream *txtout = output->txtoutput ;
		
	*txtout << "   -- integrated true infos bank  (" << thisHitBank.idtag + RAWINT_ID << ", 0) -- " << endl;
	for(map<int, string>::iterator it =  rawBank.orderedNames.begin(); it != rawBank.orderedNames.end(); it++)
	{
		int bankId   = rawBank.getVarId(it->second);
		int bankType = rawBank.getVarBankType(it->second);
		
		// we only need the first hit to get the definitions
		map<string, double> raws = HO[0].getRaws();

		// bankID 0 is hit index
		if(raws.find(it->second) != raws.end() && bankId >= 0 && bankType == RAWINT_ID)
		{
			*txtout << "    - (" << rawBank.idtag + thisHitBank.idtag << ", " << bankId << ") " << it->second << ":\t" ;
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, double> theseRaws = HO[nh].getRaws();
				*txtout <<  std::setprecision(12) << theseRaws[it->second] << "\t" ;
			}
			*txtout << endl;
		}
	}
	*txtout << "   -- end of integrated raw bank." << endl;
}




// write out true information step by step. This is common to all banks
// and not contained in the banks definitions
void txt_output ::  writeG4RawAll(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank allRawsBank = getBankFromMap("allraws", banksMap);

	initBank(output, thisHitBank);
	ofstream *txtout = output->txtoutput ;
		
	*txtout << "   -- step by step true infos bank  (" << thisHitBank.idtag + RAWSTEP_ID << ", 0) -- " << endl;
	for(map<int, string>::iterator it =  allRawsBank.orderedNames.begin(); it != allRawsBank.orderedNames.end(); it++)
	{
		int bankId   = allRawsBank.getVarId(it->second);
		int bankType = allRawsBank.getVarBankType(it->second);
		
		// we only need the first hit to get the definitions
		map<string, vector<double> > allRaws = HO[0].getAllRaws();

		// bankID 0 is hit index
		if(allRaws.find(it->second) != allRaws.end() && bankId >= 0 && bankType == RAWSTEP_ID)
		{
			*txtout << "    - (" << allRawsBank.idtag + thisHitBank.idtag << ", " << bankId << ") " << it->second << ":\t" ;
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
	*txtout << "   -- end of integrated raw bank." << endl;
}


void txt_output ::  writeG4DgtIntegrated(outputContainer* output, vector<hitOutput> HO,  string hitType, map<string, gBank> *banksMap)
{
	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank dgtBank = getDgtBankFromMap(hitType, banksMap);
	
	initBank(output, thisHitBank);
	ofstream *txtout = output->txtoutput ;
	
	*txtout << "   -- integrated digitized bank  (" << thisHitBank.idtag + DGTINT_ID << ", 0) -- " << endl;

	for(map<int, string>::iterator it =  dgtBank.orderedNames.begin(); it != dgtBank.orderedNames.end(); it++)
	{
		int bankId   = dgtBank.getVarId(it->second);
		int bankType = dgtBank.getVarBankType(it->second);
		
		// we only need the first hit to get the definitions
		map<string, double> dgts = HO[0].getDgtz();

		// bankID 0 is hit index
		if(dgts.find(it->second) != dgts.end() && bankId > 0 && bankType == DGTINT_ID)
		{
			*txtout << "    - (" << dgtBank.idtag + thisHitBank.idtag << ", " << bankId << ") " << it->second << ":\t";
			
			for(unsigned int nh=0; nh<HO.size(); nh++)
			{
				map<string, double> theseDgts = HO[nh].getDgtz();
				*txtout << theseDgts[it->second] << "\t" ;
			}
			*txtout << endl;
		}
	}
	
	*txtout << "   -- End of integrated digitized bank." << endl;
}

// index 0: hit number
// index 1: step index
// index 2: charge at electronics
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
void txt_output :: writeChargeTime(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
	gBank thisHitBank    = getBankFromMap(hitType, banksMap);
	gBank chargeTimeBank = getBankFromMap("chargeTime", banksMap);

	initBank(output, thisHitBank);
	ofstream *txtout = output->txtoutput ;

	*txtout << "   -- charge time infos (as seen by electronics) bank  (" << thisHitBank.idtag + CHARGE_TIME_ID << ", 0) -- " << endl;

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
		*txtout << "    - Hit number :\t " << thisHitN[0] << " has " << thisStep.size() << " steps." << endl;

		// identifier
		*txtout << "    - Identifier:\t " ;
		for(unsigned s=0; s<thisID.size(); s++) *txtout << (int) thisID[s] << " "  ;
		*txtout << endl ;


		// step index
		*txtout << "    - Step number:\t " ;
		for(unsigned s=0; s<thisStep.size(); s++) *txtout << (int) thisStep[s] << " "  ;
		*txtout << endl ;

		// charge at electronics
		*txtout << "    - Charge at electronics:\t " ;
		for(unsigned s=0; s<thisCharge.size(); s++) *txtout << thisCharge[s] << " "  ;
		*txtout << endl ;

		// time at electronics
		*txtout << "    - Time at electronics:\t " ;
		for(unsigned s=0; s<thisTime.size(); s++) *txtout << thisTime[s] << " "  ;
		*txtout << endl ;


	}
	*txtout << "   -- end of charge time infos (as seen by electronics) bank  (" << thisHitBank.idtag + CHARGE_TIME_ID << ", 0) -- " << endl;

}

// write fadc mode 1 (full signal shape) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void txt_output :: writeFADCMode1(outputContainer* output, vector<hitOutput> HO)
{
}

// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void txt_output :: writeFADCMode7(outputContainer* output, vector<hitOutput> HO)
{
}


void txt_output :: writeEvent(outputContainer* output)
{
	if(insideBank.size())
	{
		ofstream *txtout = output->txtoutput ;
		*txtout << " ---- End of Event  ---- " << endl;
		
	}
}

