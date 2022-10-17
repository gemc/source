// G4 headers
#include "G4MaterialPropertyVector.hh"
#include "Randomize.hh"

// gemc headers
#include "ltcc_hitprocess.h"

// C++ headers
#include <set>

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static ltccConstants initializeLTCCConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	ltccConstants ltccc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return ltccc;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	// database
	ltccc.runNo = runno;
	if(getenv ("CCDB_CONNECTION") != nullptr)
		ltccc.connection = (string) getenv("CCDB_CONNECTION");
	else
		ltccc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ltccc.connection));
	
	vector<vector<double> > data;
	// layer = left or right side
	// component = segment number
	int sector, layer, component;
	
	sprintf(ltccc.database,"/calibration/ltcc/spe:%d:%s%s",ltccc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ltccc.database);
	for(unsigned row = 0; row < data.size(); row++) {
		sector    = data[row][0] - 1;
		layer     = data[row][1] - 1;
		component = data[row][2] - 1;
		
		ltccc.speMean[sector][layer][component]  = data[row][3];
		ltccc.speSigma[sector][layer][component] = data[row][4];
		
		//		cout << " Loading ltcc sector " << sector << "  side " << layer << "  segment " << component;
		//		cout << "  spe mean: " << ltccc.speMean[sector][layer][component] << "  spe sigma: " <<  ltccc.speSigma[sector][layer][component] << endl;
	}
	
	sprintf(ltccc.database,"/calibration/ltcc/time_offsets:%d:%s%s", ltccc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ltccc.database);
	for(unsigned row = 0; row < data.size(); row++) {
		sector    = data[row][0] - 1;
		layer     = data[row][1] - 1;
		component = data[row][2] - 1;
		
		ltccc.timeOffset[sector][layer][component]  = data[row][3];
		ltccc.timeRes[sector][layer][component]  = data[row][4];
		
		//		cout << " Loading ltcc sector " << sector << "  side " << layer << "  segment " << component;
		//		cout << "  spe mean: " << ltccc.speMean[sector][layer][component] << "  spe sigma: " <<  ltccc.speSigma[sector][layer][component] << endl;
	}
	
	return ltccc;
}

map<string, double> ltcc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	int idsector  = identity[0].id;
	int idside    = identity[1].id;
	int idsegment = identity[2].id;
	int thisPid   = aHit->GetPID();
	
	if(aHit->isBackgroundHit == 1) {
		
		// background hit has all the nphe in the charge infp. Time is also first step
		double nphe     = aHit->GetCharge();
		double stepTime = aHit->GetTime()[0];
		
		dgtz["hitn"]   = hitn;
		
		dgtz["sector"]    = idsector;
		dgtz["layer"]     = idside;
		dgtz["component"] = idsegment;
		dgtz["ADC_order"] = 0;
		dgtz["ADC_ADC"]   = nphe*ltccc.speMean[idsector-1][idside-1][idsegment-1];;
		dgtz["ADC_time"]  = (int) (stepTime*24.0/1000);
		dgtz["ADC_ped"]   = 0;
		
		dgtz["TDC_order"] = 0;
		dgtz["TDC_TDC"]   = (int) (stepTime*24.0/1000);
		
		
		return dgtz;
	}
	
	
	trueInfos tInfos(aHit);
	
	// if anything else than a photon hits the PMT
	// the nphe is the particle id
	// and identifiers are negative
	// this should be changed, what if we still have a photon later?
	// if the particle is not an opticalphoton return bank filled with negative identifiers
	if(thisPid != MHit::OPTICALPHOTONPID) {
		
		dgtz["sector"]    = -idsector;
		dgtz["layer"]     = -idside;
		dgtz["component"] = -idsegment;
		
		return dgtz;
	}
	
	vector<int> tids = aHit->GetTIds();      // track ID at EACH STEP
	vector<int> pids = aHit->GetPIDs();      // particle ID at EACH STEP
	vector<double> Energies = aHit->GetEs(); // energy of the photon as it reach the pmt
	
	
	map<int, double> penergy;  // key is track id
	
	for(unsigned int s=0; s<tids.size(); s++)
	{
		// only insert the first step of each track
		// (i.e. if the map is empty
		if(penergy.find(tids[s]) == penergy.end())
			penergy[tids[s]] = Energies[s];
	}
	
	int narrived  = 0;
	int ndetected = 0;
	
	// If the detector corresponding to this hit has a material properties table with "Efficiency" defined:
	G4MaterialPropertiesTable* MPT = aHit->GetDetector().GetLogical()->GetMaterial()->GetMaterialPropertiesTable();
	G4MaterialPropertyVector* efficiency = nullptr;
	
	bool gotefficiency = false;
	if( MPT != nullptr ) {
		efficiency = (G4MaterialPropertyVector*) MPT->GetProperty("EFFICIENCY");
		if( efficiency != nullptr ) gotefficiency = true;
	}
	
	for( unsigned int iphoton = 0; iphoton<penergy.size(); iphoton++ )
	{
		//loop over all unique photons contributing to the hit:
		if( gotefficiency )
		{
			// If the material of this detector has a material properties table
			// with "EFFICIENCY" defined, then "detect" this photon with probability = efficiency
			bool outofrange = false;
			if( G4UniformRand() <= efficiency->GetValue( penergy[tids[iphoton]], outofrange ) )
				ndetected++;
			
			narrived++;
			
			if( verbosity > 4 )
			{
				cout << log_msg << " Found efficiency definition for material "
				<< aHit->GetDetector().GetLogical()->GetMaterial()->GetName()
				<< ": (Ephoton, efficiency)=(" << penergy[tids[iphoton]] << ", "
				<< ( (G4MaterialPropertyVector*) efficiency )->GetValue( penergy[tids[iphoton]], outofrange )
				<< ")" << endl;
			}
		}
		else
		{
			// No efficiency definition, "detect" all photons
			ndetected++;
		}
	}
	
	double adc = G4RandGauss::shoot(ndetected*ltccc.speMean[idsector-1][idside-1][idsegment-1], ndetected*ltccc.speSigma[idsector-1][idside-1][idsegment-1]);
	
	double timeOffset = G4RandGauss::shoot(ltccc.timeOffset[idsector-1][idside-1][idsegment-1], ltccc.timeRes[idsector-1][idside-1][idsegment-1]);
	double time = tInfos.time + timeOffset;
	
	dgtz["hitn"]   = hitn;
	
	dgtz["sector"]    = idsector;
	dgtz["layer"]     = idside;
	dgtz["component"] = idsegment;
	dgtz["ADC_order"] = 0;
	dgtz["ADC_ADC"]   = adc;
	dgtz["ADC_time"]  = time;
	dgtz["ADC_ped"]   = 0;
	
	dgtz["TDC_order"] = 0;
	dgtz["TDC_TDC"]   = (int) time;
	
	
	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}


vector<identifier>  ltcc_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}



// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ltcc_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;
	
	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);
	
	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)
	
	return noiseHits;
}


map< string, vector <int> >  ltcc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > ltcc_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double ltcc_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

void ltcc_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(ltccc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ltccc = initializeLTCCConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		ltccc.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
ltccConstants ltcc_HitProcess::ltccc = initializeLTCCConstants(-1);





