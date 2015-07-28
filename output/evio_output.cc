// gemc headers
#include "evio_output.h"
#include "utils.h"

// C++ headers
#include <fstream>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


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
	
	for(unsigned i=0; i<MAXP && i<MGP.size(); i++)
	{
		pid.push_back(MGP[i].PID);
		px.push_back(MGP[i].momentum.getX()/MeV);
		py.push_back(MGP[i].momentum.getY()/MeV);
		pz.push_back(MGP[i].momentum.getZ()/MeV);
		vx.push_back(MGP[i].vertex.getX()/mm);
		vy.push_back(MGP[i].vertex.getY()/mm);
		vz.push_back(MGP[i].vertex.getZ()/mm);
	}
	
	// creating and inserting generated particles bank  >> TAG=10 NUM=0 <<
	evioDOMNodeP generatedp = evioDOMNode::createEvioDOMNode(GENERATED_PARTICLES_BANK_TAG, 0);
	
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("pid"), bank.getVarType("pid"), pid);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("px"),  bank.getVarType("px"),  px);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("py"),  bank.getVarType("py"),  py);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("pz"),  bank.getVarType("pz"),  pz);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("vx"),  bank.getVarType("vx"),  vx);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("vy"),  bank.getVarType("vy"),  vy);
	*generatedp << addVector(GENERATED_PARTICLES_BANK_TAG, bank.getVarId("vz"),  bank.getVarType("vz"),  vz);
	
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

