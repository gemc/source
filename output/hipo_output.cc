// gemc headers
#include "hipo_output.h"
#include "utils.h"

// C++ headers
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// hipo4
#include "hipo4/writer.h"


// record the simulation conditions
// the format is a string for each variable
// the num is 0 for the mother bank, 1 for the data
void hipo_output :: recordSimConditions(outputContainer* output, map<string, string> sims)
{
	vector<string> data;
	vector<string> jdata;

//	// writing both key and argument as one string
//	for(map<string, string>::iterator it=sims.begin(); it!=sims.end(); it++) {
//		if(it->first != "JSON") {
//			data.push_back(it->first + ":  " + it->second + "  ");
//		}
//	}
}

// instantiates hipo event
// write run::config bank
void hipo_output :: writeHeader(outputContainer* output, map<string, double> data, gBank bank)
{
	if(outEvent == nullptr) {
		outEvent = new hipo::event();
	} else {
		outEvent->reset();
	}



	// Create runConfigBank with 1 row based on schema
	hipo::bank runConfigBank(output->hipoSchema->runConfigSchema, 1);

	runConfigBank.putInt("run",   0, data["runNo"]);
	runConfigBank.putInt("event", 0, data["evn"]);

	time_t t = std::time(0);
	int now = static_cast<int> (t);
	runConfigBank.putInt("unixtime", 0, now);
	runConfigBank.putInt("trigger",  0, 0);

	runConfigBank.show();

	outEvent->addStructure(runConfigBank);

}

// write user infos header
void hipo_output :: writeUserInfoseHeader(outputContainer* output, map<string, double> data)
{

}


void hipo_output :: writeRFSignal(outputContainer* output, FrequencySyncSignal rfsignals, gBank bank)
{


}

void hipo_output :: writeGenerated(outputContainer* output, vector<generatedParticle> MGP, map<string, gBank> *banksMap, vector<userInforForParticle> userInfo)
{
	double MAXP             = output->gemcOpt.optMap["NGENP"].arg;
	double SAVE_ALL_MOTHERS = output->gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg ;
	int fastMCMode          = output->gemcOpt.optMap["FASTMCMODE"].arg;  

	if(fastMCMode>0) SAVE_ALL_MOTHERS = 1;
 	if (output->gemcOpt.optMap["SAVE_ALL_ANCESTORS"].arg && (SAVE_ALL_MOTHERS == 0))
	  SAVE_ALL_MOTHERS = 1;

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

		for(unsigned int i=0; i<MAXP && i<MGP.size(); i++) {
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

	}
}

void hipo_output :: writeAncestors (outputContainer* output, vector<ancestorInfo> ainfo, gBank bank)
{
  vector<int> pid;
  vector<int> tid;
  vector<int> mtid;
  vector<double> trackE;
  vector<double> px;
  vector<double> py;
  vector<double> pz;
  vector<double> vx;
  vector<double> vy;
  vector<double> vz;
  
  for (unsigned i = 0; i < ainfo.size(); i++)
    {
      pid.push_back (ainfo[i].pid);
      tid.push_back (ainfo[i].tid);
      mtid.push_back (ainfo[i].mtid);
      trackE.push_back (ainfo[i].trackE);
      px.push_back (ainfo[i].p.getX()/MeV);
      py.push_back (ainfo[i].p.getY()/MeV);
      pz.push_back (ainfo[i].p.getZ()/MeV);
      vx.push_back (ainfo[i].vtx.getX()/MeV);
      vy.push_back (ainfo[i].vtx.getY()/MeV);
      vz.push_back (ainfo[i].vtx.getZ()/MeV);
    }


}

void hipo_output :: initBank(outputContainer* output, gBank thisHitBank, int what)
{

}


void hipo_output :: writeG4RawIntegrated(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
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
		}
	}
}


void hipo_output :: writeG4DgtIntegrated(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
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
		}
	}
}


// index 0: hit number
// index 1: step index
// index 2: charge at electronics
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
void hipo_output :: writeChargeTime(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
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


}


void hipo_output :: writeG4RawAll(outputContainer* output, vector<hitOutput> HO, string hitType, map<string, gBank> *banksMap)
{
}

void hipo_output :: writeFADCMode1( map<int,  vector<hitOutput> > HO , int ev_number)
{
}


void hipo_output :: writeFADCMode1(outputContainer* output, vector<hitOutput> HO, int ev_number)
{
}


// write fadc mode 7 (integrated mode) - jlab hybrid banks. This uses the translation table to write the crate/slot/channel
void hipo_output :: writeFADCMode7(outputContainer* output, vector<hitOutput> HO, int ev_number)
{
}

void hipo_output :: writeEvent(outputContainer* output)
{
	output->hipoWriter->addEvent(*outEvent);
}

