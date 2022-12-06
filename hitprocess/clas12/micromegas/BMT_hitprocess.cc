// gemc headers
#include "BMT_hitprocess.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Randomize.hh"
#include "G4ChargedGeantino.hh"
#include "G4Poisson.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;


static bmtConstants initializeBMTConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	bmtConstants bmtc;
	if(runno == -1) return bmtc;
	
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	
	bmtc.runNo = runno;
	
	if(getenv ("CCDB_CONNECTION") != nullptr) {
		bmtc.connection = (string) getenv("CCDB_CONNECTION");
	} else {
		bmtc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	}
	
	vector<vector<double> > data;
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(bmtc.connection));
	
	// Load the geometrical constant for each layer
	sprintf(bmtc.database,"/geometry/cvt/mvt/bmt_layer_noshim:%d:%s%s", bmtc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,bmtc.database);
	
	for(unsigned row = 0; row < data.size(); row++) {
		bmtc.AXIS[row] = data[row][1];
		bmtc.RADIUS[row] = data[row][3];
		bmtc.ZMIN[row]   = data[row][4];
		bmtc.ZMAX[row] = data[row][5];
		bmtc.EDGE1[row] = (data[row][6]-30)*degree;
		bmtc.EDGE2[row] = (data[row][7]-30)*degree;
		bmtc.NSTRIPS[row] = data[row][8];
		bmtc.hDrift =  data[row][9];
	}
	
	// Load the strip structure of each layer, compute PITCH
	bmtc.GROUP.resize(bmtc.NLAYERS); // 6 Layers
	bmtc.PITCH.resize(bmtc.NLAYERS);
	
	for (int layer=0; layer<bmtc.NLAYERS;layer++){
		sprintf(bmtc.database,"/geometry/cvt/mvt/bmt_strip_L%d:%d:%s%s", layer+1, bmtc.runNo, digiVariation.c_str(), timestamp.c_str());
		data.clear(); calib->GetCalib(data,bmtc.database);
		
		bmtc.GROUP[layer].resize(data.size());
		bmtc.PITCH[layer].resize(data.size());
		for(unsigned row = 0; row < data.size(); row++)
		{
			bmtc.GROUP[layer][row] = data[row][0];
			bmtc.PITCH[layer][row] = data[row][1];
			
			if (bmtc.AXIS[layer]==1) {//Compute angular pitch and redefine the phi angular coverage to be consistent with the pitch
				bmtc.PITCH[layer][row] = bmtc.PITCH[layer][row]/bmtc.RADIUS[layer]; //Get an angular pitch
				for (int j = 0; j <bmtc.NSECTORS ; ++j)
				{
					double middle=(bmtc.EDGE1[layer]+bmtc.EDGE2[layer])/2.;
					bmtc.EDGE1[layer] = middle-bmtc.GROUP[layer][row]*bmtc.PITCH[layer][row]/2.;
					bmtc.EDGE2[layer] = middle+bmtc.GROUP[layer][row]*bmtc.PITCH[layer][row]/2.;
				}
			}
		}
		
		for (int j = 0; j <bmtc.NSECTORS ; ++j)
		{
			if (bmtc.AXIS[layer]==1) bmtc.HV_DRIFT[layer][j]=1800;
			if (bmtc.AXIS[layer]==0) bmtc.HV_DRIFT[layer][j]=1500;
			bmtc.HV_STRIPS[layer][j]=520;
		}
	}
	
	// all dimensions are in mm
	bmtc.SigmaDrift = 0.036; //mm-1
	bmtc.hStrip2Det = bmtc.hDrift/2.;
	bmtc.nb_sigma=4;
	//bmtc.changeFieldScale(-1);  // this needs to be read from DB
	
	bmtc.Lor_Angle.Initialize(runno);
	
	// get hit time distribution parameters
	sprintf(bmtc.database,"/calibration/mvt/bmt_time:%d:%s%s", bmtc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,bmtc.database);
	bmtc.Twindow  = data[0][3]*ns;
	bmtc.Tmean    = data[0][4]*ns;
	bmtc.Tsigma   = data[0][5]*ns;
	
	// now connecting to target geometry to get its position
	sprintf(bmtc.database,"/geometry/target:%d:%s%s", bmtc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,bmtc.database);
	bmtc.targetZPos = data[0][3]*cm;
	
	return bmtc;
}

map<string, double>  BMT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double>  dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	if(aHit->isBackgroundHit == 1) {
		
		// background hit has all the energy in the first step.
		double totEdep = aHit->GetEdep()[0];
		
		int sector = identity[0].id;
		int layer  = identity[1].id;
		int strip  = identity[2].id;
		
		dgtz["sector"]    = sector;
		dgtz["layer"]     = layer;
		dgtz["component"] = strip;  // strip number
		dgtz["ADC_order"] = 0;
		dgtz["ADC_ADC"]   = int(1e6*totEdep/bmtc.w_i);
		dgtz["ADC_time"]  = bmtc.Twindow*G4UniformRand();;
		dgtz["ADC_ped"]   = 0;
		
		return dgtz;
	}
	
	// BMT ID:
	// layer, type, sector, strip
	
	int layer  = identity[0].id;
	int sector = identity[2].id;
	int strip  = identity[3].id;
	trueInfos tInfos(aHit);
	
	if(verbosity>4) {
		trueInfos tInfos(aHit);
		cout <<  log_msg << " layer: " << layer << "  sector: " << sector << "  Strip: " << strip
		<< " x=" << tInfos.x << " y=" << tInfos.y << " z=" << tInfos.z << endl;
	}
	
	
	
	dgtz["sector"]    = sector;
	dgtz["layer"]     = layer;
	dgtz["component"] = strip;  // strip number
	dgtz["ADC_order"] = 0;
	dgtz["ADC_ADC"]   = int(1e6*tInfos.eTot/bmtc.w_i);
	
	// for geantinos, assigning ADC = 1
	// notice for geant4 < 10.7, the pid of 0 conflicts with the optical photon
	if (aHit->GetPID() == 0) {
		dgtz["ADC_ADC"]   = 1.0;
	}
	
	dgtz["ADC_time"]  = bmtc.Tmean+bmtc.Tsigma*G4RandGauss::shoot(0., 1.0);
	dgtz["ADC_ped"]   = 0;
	
	
	if (strip==-1) {
		dgtz["ADC_ADC"]   = 0;
	}
	
	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}



vector<identifier>  BMT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	class bmt_strip bmts;
	
	int layer  = yid[0].id;
	int sector = yid[2].id;
	G4ThreeVector   xyz  = aStep->GetPostStepPoint()->GetPosition();  ///< Global Coordinates of interaction
	G4ThreeVector  lxyz  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz); ///< Local Coordinates of interaction
	double z0 = bmtc.ZMIN[layer-1]+Detector.dimensions[2];
	lxyz.setZ(lxyz.z()+z0);
	// if the scale is not set, then use fieldmanager to get the value
	// if fieldmanager is not found, the field is zero
	/*	if(bmtc.fieldScale == -1)
	 {*/
	const double point[4] = {xyz.x(), xyz.y(), xyz.z(), 10};
	double fieldValue[3] = {0, 0, 0};
	double phi_p=atan2(xyz.y(),xyz.x());
	
	G4FieldManager *fmanager = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetFieldManager();
	
	G4ThreeVector dm_Z(-sin(phi_p),cos(phi_p),0); //Unit vector indicating the direction of the measurement
	G4ThreeVector dm_C(0,0,1); //Unit vector indicating the direction of measurement
	
	// if no field manager, the field is zero
	if(fmanager) {
		fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);
		G4ThreeVector BField(fieldValue[0],fieldValue[1],fieldValue[2]);
		G4ThreeVector qEField(cos(phi_p),sin(phi_p),0); //Product qE
		G4ThreeVector Fdir=qEField.cross(BField); //Direction of lorentz drift
		bmtc.ThetaL=bmtc.Lor_Angle.GetAngle(bmtc.HV_DRIFT[layer-1][sector-1]/bmtc.hDrift*10,BField.perp(qEField)/gauss/1000.)*degree;
		bmtc.Theta_Ls_Z=Fdir.angle(dm_Z);
		bmtc.Theta_Ls_C=dm_C.angle(Fdir);
		
		if(bmtc.runNo == 0){
			cout << " > BMT: Field found with value " << fieldValue[2]/gauss << " gauss. Setting Lorentz angle accordingly." << endl;
			bmtc.ThetaL=bmtc.Lor_Angle.GetAngle(bmtc.HV_DRIFT[layer-1][sector-1]/bmtc.hDrift*10,BField.perp(qEField)/gauss/1000.)*degree;
			bmtc.Theta_Ls_Z=Fdir.angle(dm_Z);
			bmtc.Theta_Ls_C=dm_C.angle(Fdir);
		}
	} else {
		bmtc.ThetaL=0;
		bmtc.Theta_Ls_Z=0;
		bmtc.Theta_Ls_C=0;
		if(bmtc.runNo == 0)
			cout << " > BMT: No field found. Lorentz angle set to zero." << endl;
	}
	
	double depe = aStep->GetTotalEnergyDeposit();
	
	// resetting depe for geantinos
	if (aStep->GetTrack()->GetDefinition() == G4ChargedGeantino::ChargedGeantinoDefinition() ){
		int np=G4Poisson( (aStep->GetStepLength()/cm) *1e4); // Warning... StepLength must be in cm... because it is 10 e- per cm for MIP
		int nsec=0;
		for (int n=0;n<np;n++) nsec+=G4Poisson(2); // For each primary ionization, there are 2 secondary ionization on average.
		depe=(np+nsec)*bmtc.w_i/1e6; // So that we recover np+nsec in bmt_strip.cc line 28
	}
	
	//cout << "resolMM " << layer << " " << xyz.x() << " " << xyz.y() << " " << xyz.z() << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;
	
	// pairs [strip ID / fraction of energy ]
	vector<double> multi_hit = bmts.FindStrip(layer, sector, lxyz, depe, bmtc); //return strip=-1 and signal -1 if depe<ionization
	
	int n_multi_hits = multi_hit.size()/2;
	
	
	// sector = identity[0].id;
	// layer  = identity[1].id;
	// strip  = identity[2].id;
	
	// multi_hit[0] = id
	// multi_hit[1] = sharing
	
	// setting strip id to the one returned by the first entry in FindStrip
	yid[3].id = (int) multi_hit[0];
	
	// setting geantinoDepe
	yid[0].geantinoDepe = depe;
	
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
	if (multi_hit[1]!=-1) yid[3].id_sharing = multi_hit[1]/(1.0*(int) (1e6*depe/bmtc.w_i));
	
	// adding the additional identifiers
	for(int h=1; h<n_multi_hits; h++) {
		for(int j=0; j<3; j++) {
			identifier this_id;
			this_id.name       = yid[j].name;
			this_id.rule       = yid[j].rule;
			this_id.id         = yid[j].id;
			this_id.time       = yid[j].time;
			this_id.TimeWindow = yid[j].TimeWindow;
			this_id.TrackId    = yid[j].TrackId;
			this_id.id_sharing = multi_hit[2*h+1]; //Number of electron collected by the strip for this step
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
		this_id.id_sharing = multi_hit[2*h+1]/(1.0*(int) (1e6*depe/bmtc.w_i));//Fraction of depe: Mcoll/Nel
		yid.push_back(this_id);
	}
	
	
	return yid;
}


map< string, vector <int> >  BMT_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> BMT_HitProcess :: electronicNoise()
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
map< int, vector <double> > BMT_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double BMT_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}


void BMT_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(bmtc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		bmtc = initializeBMTConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		bmtc.runNo = runno;
	}
}


// this static function will be loaded first thing by the executable
bmtConstants BMT_HitProcess::bmtc = initializeBMTConstants(1);








