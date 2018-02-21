// gemc headers
#include "BMT_hitprocess.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "CLHEP/Vector/ThreeVector.h"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;


static bmtConstants initializeBMTConstants(int runno)
{
	// all these constants should be read from CCDB
	bmtConstants bmtc;
	bmtc.runNo = runno;
	
	if(getenv ("CCDB_CONNECTION") != NULL) {
		bmtc.connection = (string) getenv("CCDB_CONNECTION");
	} else {
		bmtc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	}
	
	bmtc.variation  = "default";
	vector<vector<double> > data;
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(bmtc.connection));

	//Load the geometrical constant for each layer
	sprintf(bmtc.database,"/geometry/cvt/mvt/bmt_layer");
	data.clear(); calib->GetCalib(data,bmtc.database);

	for(unsigned row = 0; row < data.size(); row++)
	{
	  bmtc.AXIS[row] = data[row][1];
	  bmtc.RADIUS[row] = data[row][3];
	  bmtc.ZMIN[row]   = data[row][4];
	  bmtc.ZMAX[row] = data[row][5];
	  bmtc.EDGE1[row][0] = (data[row][6]+120.)*degree;
	  bmtc.EDGE2[row][0] = (data[row][7]+120.)*degree;
	  bmtc.EDGE1[row][1] = (data[row][6])*degree;
	  bmtc.EDGE2[row][1] = (data[row][7])*degree;
	  bmtc.EDGE1[row][2] = (data[row][6]+240.)*degree;
	  bmtc.EDGE2[row][2] = (data[row][7]-120.)*degree;
	  bmtc.NSTRIPS[row] = data[row][8];
	  bmtc.hDrift =  data[row][9];
	}

	//Load the strip structure of each layer, compute PITCH
	bmtc.GROUP.resize(bmtc.NLAYERS); //6 Layers
	bmtc.PITCH.resize(bmtc.NLAYERS);
	 
	for (int layer=0; layer<bmtc.NLAYERS;layer++){
	  sprintf(bmtc.database,"/geometry/cvt/mvt/bmt_strip_L%d",layer+1);
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
		    if (bmtc.EDGE2[layer][j]>bmtc.EDGE1[layer][j]) {
		      double middle=(bmtc.EDGE1[layer][j]+bmtc.EDGE2[layer][j])/2.;
		      bmtc.EDGE1[layer][j] = middle-bmtc.GROUP[layer][row]*bmtc.PITCH[layer][row]/2.;
		      bmtc.EDGE2[layer][j] = middle+bmtc.GROUP[layer][row]*bmtc.PITCH[layer][row]/2.;
		    }
		    if (bmtc.EDGE2[layer][j]<bmtc.EDGE1[layer][j]) {
		      double middle=(bmtc.EDGE1[layer][j]+bmtc.EDGE2[layer][j]+2*pi)/2.;
		      bmtc.EDGE1[layer][j] = middle-bmtc.GROUP[layer][row]*bmtc.PITCH[layer][row]/2.;
		      bmtc.EDGE2[layer][j] = middle+bmtc.GROUP[layer][row]*bmtc.PITCH[layer][row]/2.;
		      bmtc.EDGE2[layer][j] -=2*pi;
		    }
		    	    
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
	
	if(runno == -1)
	{
		cout << " > bmt pre-initizialization. " << endl;
		return bmtc;
	}

	return bmtc;
}

map<string, double>  BMT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double>  dgtz;
	vector<identifier> identity = aHit->GetId();

	if(aHit->isBackgroundHit == 1) {

		vector<double> eDep = aHit->GetEdep();

		// background hit has all the energy in the first step
		double totEdep = eDep[0];

		int sector = identity[0].id;
		int layer  = identity[1].id;
		int strip  = identity[2].id;

		dgtz["hitn"]   = hitn;
		dgtz["sector"] = sector;
		dgtz["layer"]  = layer;
		dgtz["strip"]  = strip;
		dgtz["Edep"]   = totEdep;
		dgtz["ADC"]   = int(1e6*totEdep/bmtc.w_i);

		return dgtz;
	}

	// BMT ID:
	// layer, type, sector, strip
	
	int layer  = identity[0].id;
	int sector = identity[2].id;
	int strip  = identity[3].id;
	trueInfos tInfos(aHit);

	if(verbosity>4)
	  {
	  trueInfos tInfos(aHit);
		cout <<  log_msg << " layer: " << layer << "  sector: " << sector << "  Strip: " << strip
			 << " x=" << tInfos.x << " y=" << tInfos.y << " z=" << tInfos.z << endl;
	  }
	
	dgtz["hitn"]   = hitn;
	dgtz["layer"]  = layer;
	dgtz["sector"] = sector;
	dgtz["strip"]  = strip;
	dgtz["Edep"]   = tInfos.eTot;
	dgtz["ADC"]   = int(1e6*tInfos.eTot/bmtc.w_i);

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



vector<identifier>  BMT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
  	vector<identifier> yid = id;
	class bmt_strip bmts;
	
	int layer  = yid[0].id;
	int sector = yid[2].id;
	G4ThreeVector   xyz  = aStep->GetPostStepPoint()->GetPosition();

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
	if(fmanager)
	  {
	    fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);
	    G4ThreeVector BField(fieldValue[0],fieldValue[1],fieldValue[2]);
	    G4ThreeVector qEField(cos(phi_p),sin(phi_p),0); //Product qE
	    G4ThreeVector Fdir=qEField.cross(BField); //Direction of lorentz drift
	    bmtc.ThetaL=bmtc.Lor_Angle.GetAngle(bmtc.HV_DRIFT[layer-1][sector-1]/bmtc.hDrift*10,BField.perp(qEField)/gauss/1000.)*degree;
	    bmtc.Theta_Ls_Z=Fdir.angle(dm_Z);
	    bmtc.Theta_Ls_C=dm_C.angle(Fdir);
	   	    
	    if(bmtc.runNo == 0){
	      cout << " > BMT: Field found with value " << fieldValue[2]/gauss << " gauss. Setting Lorentz angle accordingly." << endl;
	      bmtc.ThetaL=bmtc.ThetaL=bmtc.Lor_Angle.GetAngle(bmtc.HV_DRIFT[layer-1][sector-1]/bmtc.hDrift*10,BField.perp(qEField)/gauss/1000.)*degree;
	      bmtc.Theta_Ls_Z=Fdir.angle(dm_Z);
	      bmtc.Theta_Ls_C=dm_C.angle(Fdir);
	    }
	  }
	else
	  {
	    bmtc.ThetaL=0;
	    bmtc.Theta_Ls_Z=0;
	    bmtc.Theta_Ls_C=0;
	    if(bmtc.runNo == 0)
	      cout << " > BMT: No field found. Lorentz angle set to zero." << endl;
	  }
		
	double depe = aStep->GetTotalEnergyDeposit();
	
	//cout << "resolMM " << layer << " " << xyz.x() << " " << xyz.y() << " " << xyz.z() << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;

	vector<double> multi_hit = bmts.FindStrip(layer, sector, xyz, depe, bmtc); //return strip=-1 and signal -1 if depe<ionization
	
	int n_multi_hits = multi_hit.size()/2;
	
	// closest strip
	//yid[4].id = (int) multi_hit[0];
	yid[3].id = (int) multi_hit[0];
	
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
	if (multi_hit[1]!=-1) yid[3].id_sharing = multi_hit[1]/(1.0*(int) (1e6*depe/bmtc.w_i));
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


void BMT_HitProcess::initWithRunNumber(int runno)
{	
	if(bmtc.runNo != runno)
	{
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		bmtc = initializeBMTConstants(runno);
		bmtc.runNo = runno;
	}
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


// this static function will be loaded first thing by the executable
bmtConstants BMT_HitProcess::bmtc = initializeBMTConstants(-1);








