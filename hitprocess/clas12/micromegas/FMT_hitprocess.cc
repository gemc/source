// gemc headers
#include "FMT_hitprocess.h"
#include "fmt_strip.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4CachedMagneticField.hh"
#include "CLHEP/Vector/ThreeVector.h"


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static fmtConstants initializeFMTConstants(int runno)
{
	// all these constants should be read from CCDB
	fmtConstants fmtc;


	// database
	fmtc.runNo = runno;
	fmtc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		fmtc.connection = (string) getenv("CCDB_CONNECTION");
	else
		fmtc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	fmtc.variation  = "default";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(fmtc.connection));
	vector<vector<double> > data;
	//Load the geometrical constant for all layers
	sprintf(fmtc.database,"/geometry/fmt/fmt_global");
	data.clear(); calib->GetCalib(data,fmtc.database);
	// all dimensions are in mm

	fmtc.hDrift          = data[0][0];
	fmtc.pitch           = data[0][1];
	fmtc.R_min           = data[0][3];
	fmtc.N_str           = data[0][4]; //16 connectors * 64 strips = 1024 strips for each fmt
	fmtc.N_halfstr       = data[0][5]; //number of bottom strips in the central part
	fmtc.SigmaDrift = 0.01; //mm-1

	fmtc.w_i             = 25.0;

	
	//Load the geometrical constant for all layers
	sprintf(fmtc.database,"/geometry/fmt/fmt_layer");
	data.clear(); calib->GetCalib(data,fmtc.database);
	fmtc.Z0.resize(data.size());
	fmtc.alpha.resize(data.size());
	for(unsigned row = 0; row < data.size(); row++)
	    {
	      fmtc.Z0[row]=data[row][1];
	      fmtc.alpha[row]=data[row][2]*degree;
	      
	    }
	// Number of strips and pixels
	fmtc.N_sidestr = (fmtc.N_str-2*fmtc.N_halfstr)/2; //number of strips one side
	fmtc.y_central = fmtc.N_halfstr*fmtc.pitch/2.; // Y-limit of the central part
	fmtc.nb_sigma=4;

	fmtc.R_max = fmtc.pitch*(fmtc.N_halfstr+2*fmtc.N_sidestr)/2.;

	for (int i=0;i<fmtc.NLAYERS;i++){
	  fmtc.HV_DRIFT[i]=600;
	  fmtc.HV_STRIPS_IN[i]=520;
	  fmtc.HV_STRIPS_OUT[i]=520;
	}

	fmtc.Lor_Angle.Initialize(runno);
	
	return fmtc;
}

map<string, double>FMT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	if(aHit->isBackgroundHit == 1) {

		vector<double> eDep = aHit->GetEdep();

		double totEdep = 0;
		for(int e=0; e<eDep.size(); e++) {
			totEdep += eDep[e];
		}

		int sector = identity[0].id;
		int layer  = identity[1].id;
		int strip  = identity[2].id;

		dgtz["hitn"]   = hitn;
		dgtz["sector"] = sector;
		dgtz["layer"]  = layer;
		dgtz["strip"]  = strip;
		dgtz["Edep"]   = totEdep;
		dgtz["ADC"]   = int(1e6*totEdep/fmtc.w_i);

		return dgtz;
	}

	// FMT ID:
	// layer, type, sector, strip
	
	int layer  = 1*identity[0].id + identity[1].id - 1 ; // modified on 7/27/2015 to match new geometry (Frederic Georges)
	//int layer  = 2*identity[0].id + identity[1].id - 2 ;
	int sector = identity[2].id;
	int strip  = identity[3].id;
	trueInfos tInfos(aHit);
	if(verbosity>4)
	{
	  cout <<  log_msg << " layer: " << layer << "  sector: " << sector << "  Strip: " << strip
			 << " x=" << tInfos.x << " y=" << tInfos.y << " z=" << tInfos.z << endl;
	}
	
	dgtz["hitn"]   = hitn;
	dgtz["layer"]  = layer;
	dgtz["sector"] = sector;
	dgtz["strip"]  = strip;
	dgtz["Edep"]   = tInfos.eTot;
	dgtz["ADC"]   = (int) (tInfos.eTot*1e6/fmtc.w_i);
	
	if (strip==-1) {
	  dgtz["Edep"]   = 0;
	  dgtz["ADC"]   = 0;
	}
	
	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	bool rejectHitConditions = false;
	if(rejectHitConditions) {
	  writeHit = false;
	}

	return dgtz;
}



vector<identifier>  FMT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	double x, y, z;
	G4ThreeVector   xyz    = aStep->GetPostStepPoint()->GetPosition();
	x = xyz.x();
	y = xyz.y();
	z = xyz.z();
	double point[4] = {xyz.x(), xyz.y(), xyz.z(),10};
	double fieldValue[6] = {0, 0, 0, 0, 0, 0};

	vector<identifier> yid = id;
	class fmt_strip fmts;
		
	int layer  = 1*yid[0].id + yid[1].id - 1 ; // modified on 7/27/2015 to match new geometry (Frederic Georges)
	int sector = yid[2].id;
	
	G4FieldManager *fmanager = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetFieldManager();

		// if no field manager, the field is zero
	if(fmanager)
	  {
	    fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);
	   
	    G4ThreeVector BField(fieldValue[0],fieldValue[1],fieldValue[2]);
	    G4ThreeVector qEField(0,0,1); //Product q*v
	    G4ThreeVector Fdir=qEField.cross(BField); //Direction of lorentz drift
	    fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);
	    fmtc.ThetaL=fmtc.Lor_Angle.GetAngle(fmtc.HV_DRIFT[layer-1]/fmtc.hDrift*10,BField.perp(qEField)/gauss/1000.)*degree;
	    fmtc.Theta_Ls=atan2(Fdir.y(),Fdir.x());
	    	    	   	    
	    if(fmtc.runNo == 0){
	      cout << " > BMT: Field found with value " << fieldValue[2]/gauss << " gauss. Setting Lorentz angle accordingly." << endl;
	      fmtc.ThetaL=fmtc.Lor_Angle.GetAngle(fmtc.HV_DRIFT[layer-1]/fmtc.hDrift*10,BField.perp(qEField)/gauss/1000.)*degree;
	      fmtc.Theta_Ls=atan2(Fdir.y(),Fdir.x());
	    }
	  }
	else
	  {
	    fmtc.ThetaL=0;
	    fmtc.Theta_Ls=0;
	    if(fmtc.runNo == 0) cout << " > BMT: No field found. Lorentz angle set to zero." << endl;
	  }
	
	//yid[3].id = fmts.FindStrip(layer-1, sector-1, x, y, z);
	double depe = aStep->GetTotalEnergyDeposit();
	//cout << "resolMM " << layer << " " << x << " " << y << " " << z << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;
	vector<double> multi_hit = fmts.FindStrip(layer-1, sector-1, x, y, z, depe, fmtc);
	
	int n_multi_hits = multi_hit.size()/2;
	
	// closest strip
	//yid[4].id = (int) multi_hit[0];
	yid[3].id = (int) multi_hit[0];
	
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
	if (multi_hit[1]!=-1) yid[3].id_sharing = multi_hit[1]/(1.0*(int) (1e6*depe/fmtc.w_i));
	// yid[4].id_sharing = multi_hit[1];
	
	// additional strip
	for(int h=1; h<n_multi_hits; h++)
	{
		for(int j=0; j<3; j++)
		{
			identifier this_id;
			this_id.name       = yid[j].name;
			this_id.rule       = yid[j].rule;
			this_id.id         = yid[j].id;
			this_id.time       = yid[j].time;
			this_id.TimeWindow = yid[j].TimeWindow;
			this_id.TrackId    = yid[j].TrackId;
			this_id.id_sharing = multi_hit[2*h+1];
			yid.push_back(this_id);
		}
		// last id is strip
		identifier this_id;
		this_id.name       = yid[3].name;
		this_id.rule       = yid[3].rule;
		this_id.id         = (int) multi_hit[2*h];
		this_id.time       = yid[3].time;
		this_id.TimeWindow = yid[3].TimeWindow;
		this_id.TrackId    = yid[3].TrackId;
		this_id.id_sharing = multi_hit[2*h+1]/(1.0*(int) (1e6*depe/fmtc.w_i));
		yid.push_back(this_id);
	}
	
	return yid;
}

void FMT_HitProcess::initWithRunNumber(int runno)
{
	if(fmtc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		fmtc = initializeFMTConstants(runno);
		fmtc.runNo = runno;
	}
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> FMT_HitProcess :: electronicNoise()
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
map< int, vector <double> > FMT_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double FMT_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

map< string, vector <int> >  FMT_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


// this static function will be loaded first thing by the executable
fmtConstants FMT_HitProcess::fmtc = initializeFMTConstants(-1);











