// G4 headers
#include "G4MaterialPropertyVector.hh"
#include "Randomize.hh"

// gemc headers
#include "htcc_hitprocess.h"
#include "detector.h"
#include <set>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static htccConstants initializeHTCCConstants(int runno)
{
	// all these constants should be read from CCDB
	htccConstants htccc;

	// do not initialize at the beginning, only after the end of the first event,                                                                                       
        // with the proper run number coming from options or run table                                                                                                      
        if(runno == -1) return htccc;

	// database
	htccc.runNo = runno;
	htccc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		htccc.connection = (string) getenv("CCDB_CONNECTION");
	else
		htccc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	htccc.variation  = "main";
	int isec,ilay,istr;

        vector<vector<double> > data;

	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(htccc.connection));
	cout<<"HTCC:Getting status"<<endl;
        sprintf(htccc.database,"/calibration/htcc/status:%d",htccc.runNo);
        data.clear() ; calib->GetCalib(data,htccc.database);
        for(unsigned row = 0; row < data.size(); row++)
	  {
	    isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	    htccc.status[isec-1][ilay-1].push_back(data[row][3]);
	  }

	cout<<"HTCC:Getting time_offset"<<endl;
        sprintf(htccc.database,"/calibration/htcc/time_offset:%d",htccc.runNo);
        data.clear() ; calib->GetCalib(data,htccc.database);
        for(unsigned row = 0; row < data.size(); row++)
	  {
	    isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	    htccc.tshift[isec-1][ilay-1].push_back(data[row][3]);

	  }

	return htccc;
}

map<string, double> htcc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	
    // we want to crash if identity doesn't have size 3
    vector<identifier> identity = aHit->GetId();
    int idsector = identity[0].id;
    int idring   = identity[1].id;
    int idhalf   = identity[2].id;
    int thisPid  = aHit->GetPID();
    
    trueInfos tInfos(aHit);
    
    int ndetected;
    // if anything else than a photon hits the PMT
    // the nphe is the particle id
    // and identifiers are negative
	 // this should be changed, what if we still have a photon later?
    dgtz["sector"] = -idsector;
    dgtz["ring"]   = -idring;
    dgtz["half"]   = -idhalf;
    dgtz["nphe"]   = thisPid;
    dgtz["time"]   = tInfos.time;
    dgtz["hitn"]   = hitn;

    
	// if the particle is not an opticalphoton return bank filled with negative identifiers
	if(thisPid != 0)
        return dgtz;
	
	// Since the HTCC hit involves a PMT which detects photons with a certain quantum efficiency (QE)
	// we want to implement QE here in a flexible way:
	
	// That means we want to find out if the material of the photocathode has been defined,
	// and if so, find out if it has a material properties table which defines an efficiency.
	// if we find both of these properties, then we accept this event with probability QE:
	
	vector<int> tids = aHit->GetTIds(); // track ID at EACH STEP
	vector<int> pids = aHit->GetPIDs(); // particle ID at EACH STEP
	set<int> TIDS;                      // an array containing all UNIQUE tracks in this hit
	vector<double> photon_energies;
	
	vector<double> Energies = aHit->GetEs();
	
	
	// this needs to be optimized
	// uaing the return value of insert is unnecessary

	for(unsigned int s=0; s<tids.size(); s++)
	{
		// selecting optical photons
		if(pids[s] == 0)
		{
			// insert this step into the set of track ids (set can only have unique values).
			pair< set<int> ::iterator, bool> newtrack = TIDS.insert(tids[s]);
			
			// if we found a new track, then add the initial photon energy of this
			// track to the list of photon energies, for when we calculate quantum efficiency later
			if( newtrack.second ) photon_energies.push_back( Energies[s] );
		}
	}
	
	
	// here is the fun part: figure out the number of photons we detect based
	// on the quantum efficiency of the photocathode material, if defined:
	
	
	// If the detector corresponding to this hit has a material properties table with "Efficiency" defined:
	G4MaterialPropertiesTable* MPT = aHit->GetDetector().GetLogical()->GetMaterial()->GetMaterialPropertiesTable();
	G4MaterialPropertyVector* efficiency = NULL;
	ndetected = 0;
	bool gotefficiency = false;
	if( MPT != NULL )
	{
		efficiency = (G4MaterialPropertyVector*) MPT->GetProperty("EFFICIENCY");
		if( efficiency != NULL ) gotefficiency = true;
	}
	
	for( unsigned int iphoton = 0; iphoton<TIDS.size(); iphoton++ )
	{
		//loop over all unique photons contributing to the hit:
		if( gotefficiency )
		{
			// If the material of this detector has a material properties table
			// with "EFFICIENCY" defined, then "detect" this photon with probability = efficiency
			bool outofrange = false;
			if( G4UniformRand() <= efficiency->GetValue( photon_energies[iphoton], outofrange ) )
				ndetected++;
			
			if( verbosity > 4 )
			{
				cout << log_msg << " Found efficiency definition for material "
				<< aHit->GetDetector().GetLogical()->GetMaterial()->GetName()
				<< ": (Ephoton, efficiency)=(" << photon_energies[iphoton] << ", "
				<< ( (G4MaterialPropertyVector*) efficiency )->GetValue( photon_energies[iphoton], outofrange )
				<< ")" << endl;
			}
		}
		else
		{
			// No efficiency definition, "detect" all photons
			ndetected++;
		}
	}
	
	if(verbosity>4)
	{
	
		cout <<  log_msg << " (sector, ring, half)=(" << idsector << ", " << idring << ", " << idhalf << ")"
		     << " x=" << tInfos.x/cm << " y=" << tInfos.y/cm << " z=" << tInfos.z/cm << endl;
	}
	//status flags
	switch (htccc.status[idsector-1][idhalf-1][idring-1])
          {
          case 0:
            break;
          case 1:
            ndetected=0;
            break;
          case 2:
            tInfos.time = 0;

            break;
          case 3:
            tInfos.time = 0;
            ndetected = 0;

            break;
          case 5:
            break;

          default:
            cout << " > Unknown HTCC status: " << htccc.status[idsector-1][idhalf-1][idring-1] << " for sector " << idsector << ",  halfsector " << idhalf << ", ring  " << idring << endl;
          }





	dgtz["sector"] = idsector;
	dgtz["ring"]   = idring;
	dgtz["half"]   = idhalf;
	dgtz["nphe"]   = ndetected;
	dgtz["time"]   = tInfos.time + htccc.tshift[idsector-1][idhalf-1][idring-1];
	dgtz["hitn"]   = hitn;

	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	bool rejectHitConditions = false;
	if(rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
}


vector<identifier>  htcc_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

void htcc_HitProcess::initWithRunNumber(int runno)
{
	if(htccc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		htccc = initializeHTCCConstants(runno);
		htccc.runNo = runno;
	}
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> htcc_HitProcess :: electronicNoise()
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


map< string, vector <int> >  htcc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > htcc_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double htcc_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

// this static function will be loaded first thing by the executable
htccConstants htcc_HitProcess::htccc = initializeHTCCConstants(-1);







