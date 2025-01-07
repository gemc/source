// gemc headers
#include "recoil_hitprocess.h"
#include "recoil_strip.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4CachedMagneticField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "G4Trap.hh"
#include "G4Box.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static recoilConstants initializerecoilConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	recoilConstants recoilC;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return recoilC;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	// database
	recoilC.runNo = runno;
	recoilC.date       = "2024-11-29";
	/*	if(getenv ("CCDB_CONNECTION") != nullptr)
		recoil.connection = (string) getenv("CCDB_CONNECTION");
	else
	recoil.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";*/
	
	/*
	 unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(fmtc.connection));
	 vector<vector<double> > data;
	 //Load the geometrical constant for all layers
	 sprintf(fmtc.database,"/geometry/recoil/recoil_global:%d:%s%s", fmtc.runNo, digiVariation.c_str(), timestamp.c_str());
	 data.clear(); calib->GetCalib(data,recoil.database);
	 // all dimensions are in mm
	 */
	/*number of strip in each chambers*/
	/*	recoilC.number_strip_chamber_v[0] = 196;
	recoilC.number_strip_chamber_v[1] = 325;
	recoilC.number_strip_chamber_v[2] = 465;
	recoilC.number_strip_chamber_u[0] = 156;
	recoilC.number_strip_chamber_u[1] = 259;
	recoilC.number_strip_chamber_u[2] = 371;

	recoilC.number_strip_chamber[0] = 352;
	recoilC.number_strip_chamber[1] = 584;
	recoilC.number_strip_chamber[2] = 836;
	
	recoilC.number_of_strip = 1772; //Total number of strip */
	recoilC.stripU_stereo_angle = 0 ; // angle between strip and trapezoid base in degree
	recoilC.stripU_pitch = 1.;  //mm

	recoilC.stripV_width[0] = 0.4;
	recoilC.stripV_width[1] = 0.4;
	recoilC.stripV_width[2] = 0.4;
	recoilC.stripV_width[3] = 0.4;

	recoilC.stripV_stereo_angle = 90 ; // angle between strip and trapezoid base in degree
	recoilC.stripV_pitch = 1;  //mm

	recoilC.stripU_width[0] = 0.4;
	recoilC.stripU_width[1] = 0.4;
	recoilC.stripU_width[2] = 0.4;
	recoilC.stripU_width[3] = 0.4;
	
	recoilC.w_i=25; //ionization potential assumed to be 25 eV
	recoilC.sigma_td= 0.5;         // effective value to take into account transverse diffusion + charge dispersion
	recoilC.nb_sigma = 5;            // Number of sigma to study around the closest strip
	recoilC.gain =1E4;
	
	// drift velocity
	recoilC.v_drift = 5E-3; // velocity drift [cm/ns]
	recoilC.sigma_time = 20; // time resolution 20 ns
	
	return recoilC;
}

map<string, double>recoil_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	//recoilConstants recoilC;
	trueInfos tInfos(aHit);
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = identity[1].id;
	dgtz["layer"]  = identity[3].id;
	dgtz["component"]  = identity[4].id;
    if(identity[4].id ==-15000){
    	dgtz["ADC_ADC"]  = 0;
    	dgtz["ADC_time"] = 0;
    }else{
    	dgtz["ADC_ADC"]  = (1.0*(int) (recoilC.gain*1e6*tInfos.eTot/recoilC.w_i));
    	dgtz["ADC_time"] = identity[4].time;
    }
	dgtz["ADC_ped"]   = 0;
	
	//	cout<<"Sector= "<<dgtz["sector"]<<" Layer= "<<dgtz["layer"]<<" Component= "<<dgtz["component"]<<" ADC= "<<dgtz["ADC_ADC"]<<endl;
	
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
	
}

vector<identifier> recoil_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	
	//recoilConstants recoilC;
	
	vector<identifier> yid;
	recoil_strip recoil_strip;
	G4ThreeVector   xyz    = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector  lxyz    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz); ///< Local Coordinates of interaction
		
	G4VTouchable* TH = (G4VTouchable*) aStep->GetPreStepPoint()->GetTouchable();
	G4Box *Box = dynamic_cast<G4Box*>(TH->GetSolid());
	recoilC.Xhalf = Box->GetXHalfLength();
	recoilC.Yhalf = Box->GetYHalfLength();
	recoilC.Zhalf = Box->GetZHalfLength();
	
	double depe = aStep->GetTotalEnergyDeposit();
	double time = aStep->GetPostStepPoint()->GetGlobalTime();
	
	/* STRIP U */
    recoilC.get_strip_info("strip_u");
    vector<recoil_strip_found> multi_hit_u = recoil_strip.FindStrip(lxyz, depe, recoilC, time);
    int n_multi_hits_u = multi_hit_u.size();

    for(int h=0; h<n_multi_hits_u; h++){
      //  cout<<"in the loop for strip_u"<<endl;
      for(int j=0; j<5; j++)
	{
	  // j=0 region ; j=1 sector; j2 chamber; j3 layer; j4 component
	  //  cout<<"strip_u multi_hit_u dep energy "<<multi_hit_u.size()<<endl;
	  identifier this_id;
	  this_id.name       = id[j].name;
	  this_id.rule       = id[j].rule;
	  if(j==0) this_id.id = id[j].id; 
	  if(j==1) this_id.id = id[j].id;
	  if(j==2) this_id.id = id[j].id;
	  if(j==3) {
	    this_id.id = 2*id[0].id-1; 
	  }
	  this_id.time       = id[j].time;
	  
	  if(j==4){    //J==4 strip ID
	    //	    if(id[2].id>0) {
	      if(multi_hit_u.at(h).numberID ==-15000){
		this_id.id  = multi_hit_u.at(h).numberID ;
	      } else {
		//		this_id.id  = multi_hit_u.at(h).numberID + std::accumulate(recoilC.number_strip_chamber_u,recoilC.number_strip_chamber_u +id[2].id-1,0);
		this_id.id  = multi_hit_u.at(h).numberID ;
	      }
	      //}else this_id.id  = multi_hit_u.at(h).numberID;
	    this_id.time       = multi_hit_u.at(h).time;
	    //	    cout<<"ID strip u "<<this_id.id<<endl;
	  }
	  
	  this_id.TimeWindow = id[j].TimeWindow;
	  this_id.TrackId    = id[j].TrackId;
	  this_id.id_sharing = multi_hit_u.at(h).weight/(recoilC.gain*1e6*depe/recoilC.w_i);
	  
	  yid.push_back(this_id);
	}
    }
    /* STRIP V */
    recoilC.get_strip_info("strip_v");
    vector<recoil_strip_found> multi_hit_v = recoil_strip.FindStrip(lxyz, depe, recoilC, time);
    int n_multi_hits_v = multi_hit_v.size();
    
    for(int h=0; h<n_multi_hits_v; h++){
      //  cout<<"in the loop for strip_v"<<endl;
      for(int j=0; j<5; j++)
	{
	  identifier this_id;
	  this_id.name       = id[j].name;
	  this_id.rule       = id[j].rule;
	  //cout<<"strip v - j "<<j<<" id[j] "<<id[j].id<<endl;
	  if(j==0) this_id.id = id[j].id;
	  if(j==1) this_id.id = id[j].id;
	  if(j==2) this_id.id = id[j].id;
	  if(j==3) this_id.id = 2*id[0].id;
	  this_id.time       = id[j].time;
	  
	  if(j==4){    //J==4 strip ID
	    //if(id[2].id>0) {
	      if(multi_hit_v.at(h).numberID ==-15000){
		this_id.id  = multi_hit_v.at(h).numberID ;
	      } else {
		//		this_id.id  = multi_hit_v.at(h).numberID + std::accumulate(recoilC.number_strip_chamber_v,recoilC.number_strip_chamber_v +id[2].id-1,0);
		this_id.id  = multi_hit_v.at(h).numberID;
	      }
	      
	      //	    }else this_id.id  = multi_hit_v.at(h).numberID;
	    this_id.time       = multi_hit_v.at(h).time;
	    //	    cout<<"ID strip v "<<this_id.id<<endl;
	  }
	  this_id.TimeWindow = id[j].TimeWindow;
	  this_id.TrackId    = id[j].TrackId;
	  this_id.id_sharing = multi_hit_v.at(h).weight/(recoilC.gain*1e6*depe/recoilC.w_i);
	  yid.push_back(this_id);
	}
    }
    return yid;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> recoil_HitProcess :: electronicNoise()
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

// - charge: returns charge/time digitized information / step
map< int, vector <double> > recoil_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double recoil_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

map< string, vector <int> >  recoil_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

void recoil_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(recoilC.runNo != runno) {
	  //	cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		recoilC = initializerecoilConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		recoilC.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
recoilConstants recoil_HitProcess::recoilC = initializerecoilConstants(1);











