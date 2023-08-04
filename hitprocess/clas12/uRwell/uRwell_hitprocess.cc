// gemc headers
#include "uRwell_hitprocess.h"
#include "uRwell_strip.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4CachedMagneticField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "G4Trap.hh"


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static uRwellConstants initializeuRwellConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	uRwellConstants urwellC;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return urwellC;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	// database
	urwellC.runNo = runno;
	urwellC.date       = "2022-08-23";
	if(getenv ("CCDB_CONNECTION") != nullptr)
		urwellC.connection = (string) getenv("CCDB_CONNECTION");
	else
		urwellC.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	/*
	 unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(fmtc.connection));
	 vector<vector<double> > data;
	 //Load the geometrical constant for all layers
	 sprintf(fmtc.database,"/geometry/uRwell/uRwell_global:%d:%s%s", fmtc.runNo, digiVariation.c_str(), timestamp.c_str());
	 data.clear(); calib->GetCalib(data,uRwell.database);
	 // all dimensions are in mm
	 */
	/*number of strip in each chambers*/
	urwellC.number_strip_chamber[0] = 542;
	urwellC.number_strip_chamber[1] = 628;
	urwellC.number_strip_chamber[2] = 714;
	
	urwellC.number_of_strip = 1884; //Total number of strip
	urwellC.stripU_stereo_angle = -10 ; // angle between strip and trapezoid base in degree
	urwellC.stripU_pitch = 1.;  //mm

	urwellC.stripU_width[0] = 0.4;
	urwellC.stripU_width[1] = 0.4;
	urwellC.stripU_width[2] = 0.4;
	urwellC.stripU_width[3] = 0.4;

	urwellC.stripU_width_proto[0] = 0.175;
	urwellC.stripU_width_proto[1] = 0.350;
	urwellC.stripU_width_proto[2] = 0.262;
	urwellC.stripU_width_proto[3] = 0.350;

	urwellC.stripV_stereo_angle = 10 ; // angle between strip and trapezoid base in degree
	urwellC.stripV_pitch = 1.;  //mm

	urwellC.stripV_width[0] = 0.4;
	urwellC.stripV_width[1] = 0.4;
	urwellC.stripV_width[2] = 0.4;
	urwellC.stripV_width[3] = 0.4;

	urwellC.stripV_width_proto[0] = 0.355;
	urwellC.stripV_width_proto[1] = 0.650;
	urwellC.stripV_width_proto[2] = 0.5;
	urwellC.stripV_width_proto[3] = 0.650;
	
	urwellC.w_i=25; //ionization potential assumed to be 25 eV
	urwellC.sigma_td= 0.5;         // effective value to take into account transverse diffusion + charge dispersion
	urwellC.nb_sigma = 5;            // Number of sigma to study around the closest strip
	urwellC.gain =1E4;
	
	// drift velocity
	urwellC.v_drift = 5E-3; // velocity drift [cm/ns]
	urwellC.sigma_time = 20; // time resolution 20 ns
	
	return urwellC;
}




map<string, double>uRwell_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	//uRwellConstants uRwellC;
	trueInfos tInfos(aHit);
	
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = identity[1].id;
	dgtz["layer"]  = identity[3].id;
	dgtz["component"]  = identity[4].id;
    if(identity[4].id ==-15000){
    	dgtz["ADC"]  = 0;
    	dgtz["time"]   = 0;
    }else{
    	dgtz["ADC"]  = (1.0*(int) (uRwellC.gain*1e6*tInfos.eTot/uRwellC.w_i));
    	dgtz["time"]   = identity[4].time;
    }
	dgtz["ADC_ped"]   = 0;
	
//	cout<<dgtz["sector"]<<" "<<dgtz["layer"]<<" "<<dgtz["component"]<<" "<<dgtz["ADC"]<<endl;
	
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
	
}



vector<identifier> uRwell_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	
	//uRwellConstants uRwellC;
	
	vector<identifier> yid;
	uRwell_strip URwell_strip;
	// double Lorentz_angle=0;
	G4ThreeVector   xyz    = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector  lxyz    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz); ///< Local Coordinates of interaction
	
	
	
	
	G4VTouchable* TH = (G4VTouchable*) aStep->GetPreStepPoint()->GetTouchable();
	G4Trap *Trap = dynamic_cast<G4Trap*>(TH->GetSolid());
	
	//cout << Trap->GetName()<<endl;

	 bool isProto = false;
	 if(Trap->GetName().find("proto")!=std::string::npos) isProto = true;

	uRwellC.Xhalf_base = Trap->GetXHalfLength1();
	uRwellC.Xhalf_Largebase = Trap->GetXHalfLength2();
	uRwellC.Yhalf = Trap->GetYHalfLength1();
	uRwellC.Zhalf = Trap->GetZHalfLength();
	
	//int sector = id[0].id;
	// int chamber = id[1].id;
	
	double depe = aStep->GetTotalEnergyDeposit();
	double time = aStep->GetPostStepPoint()->GetGlobalTime();
	/*
	 
	 double point[4] = {xyz.x(), xyz.y(), xyz.z(),10};
	 double fieldValue[6] = {0, 0, 0, 0, 0, 0};
	 
	 G4FieldManager *fmanager = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetFieldManager();
	 
	 // if no field manager, the field is zero
	 if(fmanager)
	 {
	 fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);
	 G4ThreeVector BField(fieldValue[0],fieldValue[1],fieldValue[2]);
	 G4ThreeVector qEField(xyz_norm3.x(), xyz_norm3.y(), xyz_norm3.z()); //Product q*v
	 G4ThreeVector Fdir=qEField.cross(BField); //Direction of lorentz drift
	 cout <<"angle"<< qEField.angle(BField)*180/3.1415<<endl;
	 cout << "B field :"<< fieldValue[0]<< " "<<fieldValue[1]<< " "<<fieldValue[3]<<endl;
	 cout << "B field :"<< sqrt(fieldValue[0]*fieldValue[0] + fieldValue[1]*fieldValue[1] + fieldValue[2]*fieldValue[2])<<endl;
	 cout << "B field :"<< sqrt(fieldValue[0]*fieldValue[0] + fieldValue[1]*fieldValue[1] + fieldValue[2]*fieldValue[2])/gauss<<endl;
	 Lorentz_angle =0;
	 }
	 */
	
	/* STRIP U */
    uRwellC.get_strip_info("strip_u", isProto);
    vector<uRwell_strip_found> multi_hit_u = URwell_strip.FindStrip(lxyz, depe, uRwellC, time, isProto);
    int n_multi_hits_u = multi_hit_u.size();


    for(int h=0; h<n_multi_hits_u; h++){

    		for(int j=0; j<5; j++)
    		{
    			// j=0 region ; j=1 sector; j2 chamber; j3 layer; j4 component

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
    				if(id[2].id>0) {
					    if(multi_hit_u.at(h).numberID ==-15000){
					    	this_id.id  = multi_hit_u.at(h).numberID ;
					    } else {
					    	this_id.id  = multi_hit_u.at(h).numberID + std::accumulate(uRwellC.number_strip_chamber,uRwellC.number_strip_chamber +id[2].id-1,0);
					    }
    				}else this_id.id  = multi_hit_u.at(h).numberID;
    					this_id.time       = multi_hit_u.at(h).time;

    					}

    			this_id.TimeWindow = id[j].TimeWindow;
    			this_id.TrackId    = id[j].TrackId;
    			this_id.id_sharing = multi_hit_u.at(h).weight/(uRwellC.gain*1e6*depe/uRwellC.w_i);

    			yid.push_back(this_id);
    		}

    }

	
    /* STRIP V */
    uRwellC.get_strip_info("strip_v", isProto);
	vector<uRwell_strip_found> multi_hit_v = URwell_strip.FindStrip(lxyz, depe, uRwellC, time, isProto);
	int n_multi_hits_v = multi_hit_v.size();

	 for(int h=0; h<n_multi_hits_v; h++){
	 		for(int j=0; j<5; j++)
	 		{
	 			identifier this_id;
				this_id.name       = id[j].name;
				this_id.rule       = id[j].rule;
				if(j==0) this_id.id = id[j].id;
				if(j==1) this_id.id = id[j].id;
				if(j==2) this_id.id = id[j].id;
				if(j==3) this_id.id = 2*id[0].id;
				this_id.time       = id[j].time;

				if(j==4){    //J==4 strip ID
				   if(id[2].id>0) {
					    if(multi_hit_v.at(h).numberID ==-15000){
					    	this_id.id  = multi_hit_v.at(h).numberID ;
					    } else {
					    	this_id.id  = multi_hit_v.at(h).numberID + std::accumulate(uRwellC.number_strip_chamber,uRwellC.number_strip_chamber +id[2].id-1,0);
					    }

				   }else this_id.id  = multi_hit_v.at(h).numberID;
				   	   this_id.time       = multi_hit_v.at(h).time;
				}
				this_id.TimeWindow = id[j].TimeWindow;
				this_id.TrackId    = id[j].TrackId;
				this_id.id_sharing = multi_hit_v.at(h).weight/(uRwellC.gain*1e6*depe/uRwellC.w_i);
				yid.push_back(this_id);
	 	}
	 		}
	
	
	return yid;
	
}




// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> uRwell_HitProcess :: electronicNoise()
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
map< int, vector <double> > uRwell_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double uRwell_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

map< string, vector <int> >  uRwell_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

void uRwell_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(uRwellC.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		uRwellC = initializeuRwellConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		uRwellC.runNo = runno;
	}
}



// this static function will be loaded first thing by the executable
uRwellConstants uRwell_HitProcess::uRwellC = initializeuRwellConstants(1);











