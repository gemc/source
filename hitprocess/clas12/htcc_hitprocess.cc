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

map<string, double> htcc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// if the particle is not an opticalphoton return empty bank
	if(aHit->GetPID() != 0) return dgtz;
	
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
	unsigned nsteps = tids.size();
	
	for(unsigned int s=0; s<nsteps; s++)
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
	
	// we want to crash if identity doesn't have size 3
	int idsector = identity[0].id;
	int idring   = identity[1].id;
	int idhalf   = identity[2].id;

	
	// if identifiers for theta and phi indices have been defined, use them.
	// otherwise, calculate itheta and iphi from the (mandatory) pmt index
	// assuming the default configuration:
	int isector=-1, iring=-1, ihalf=-1;
	
	
	// what is going on here i dont think this is right
	// double check!!
	if( idsector >= 0 && idring >= 0 && idhalf >= 0 )
	{
		isector = identity[idsector].id;
		iring   = identity[idring].id;
		ihalf  = identity[idhalf].id;
	}
	else return dgtz; //no valid ID
	
	// here is the fun part: figure out the number of photons we detect based
	// on the quantum efficiency of the photocathode material, if defined:
	
	int ndetected = 0;
	
	// If the detector corresponding to this hit has a material properties table with "Efficiency" defined:
	G4MaterialPropertiesTable* MPT = aHit->GetDetector().GetLogical()->GetMaterial()->GetMaterialPropertiesTable();
	G4MaterialPropertyVector* efficiency = NULL;
	
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
			//If the material of this detector has a material properties table
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
			//No efficiency definition, "detect" all photons
			ndetected++;
		}
	}
	
	if(verbosity>4)
	{
		trueInfos tInfos(aHit);

	
		cout <<  log_msg << " (sector, ring, half)=(" << isector << ", " << iring << ", " << ihalf << ")"
		     << " x=" << tInfos.x/cm << " y=" << tInfos.y/cm << " z=" << tInfos.z/cm << endl;
	}
	dgtz["sector"] = idsector;
	dgtz["ring"]   = idring;
	dgtz["half"]   = idhalf;
	dgtz["nphe"]   = ndetected;
	dgtz["hitn"]   = hitn;
	
	return dgtz;
}


vector<identifier>  htcc_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}



map< string, vector <int> >  htcc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}








