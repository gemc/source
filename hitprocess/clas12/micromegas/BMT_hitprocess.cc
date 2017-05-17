// gemc headers
#include "BMT_hitprocess.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"

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
	}
	
	// all dimensions are in mm
	bmtc.SigmaDrift = 0.4;
	bmtc.hStrip2Det = bmtc.hDrift/2.;
	bmtc.nb_sigma=4;
	bmtc.changeFieldScale(-1);  // this needs to be read from DB

	if(runno == -1)
	{
		cout << " > bmt pre-initizialization. " << endl;
		return bmtc;
	}

	//bmtc.ZMIN[0] -= CR4C_group[j]*CR4C_width[j];
	
	return bmtc;
}

map<string, double>  BMT_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double>  dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// BMT ID:
	// layer, type, sector, strip
	
	//int layer  = 2*identity[0].id + identity[1].id - 2 ;
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
	dgtz["ADC"]   = (int) tInfos.eTot;
	return dgtz;
}



vector<identifier>  BMT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	G4ThreeVector   xyz  = aStep->GetPostStepPoint()->GetPosition();

	  // if the scale is not set, then use fieldmanager to get the value
	  // if fieldmanager is not found, the field is zero
	if(bmtc.fieldScale == -1)
	{
		const double point[4] = {xyz.x(), xyz.y(), xyz.z(), 10};
		double fieldValue[3] = {0, 0, 0};
		
		
		G4FieldManager *fmanager = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetFieldManager();
		
		// if no field manager, the field is zero
		if(fmanager)
		{
			fmanager->GetDetectorField()->GetFieldValue(point, fieldValue);
			if(bmtc.runNo == 0)
				cout << " > BMT: Field found with value " << fieldValue[2]/gauss << " gauss. Setting Lorentz angle accordingly." << endl;
			bmtc.changeFieldScale((fieldValue[2]/gauss)/50000.0);
		}
		else
		{
			bmtc.changeFieldScale(0);
			if(bmtc.runNo == 0)
				cout << " > BMT: No field found. Lorentz angle set to zero." << endl;
		}
	}
	
	vector<identifier> yid = id;
	class bmt_strip bmts;
	
	//int layer  = 2*yid[0].id + yid[1].id - 2 ;
	int layer  = yid[0].id;
	int sector = yid[2].id;
	
	double depe = aStep->GetTotalEnergyDeposit();
	//cout << "resolMM " << layer << " " << xyz.x() << " " << xyz.y() << " " << xyz.z() << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;

	vector<double> multi_hit = bmts.FindStrip(layer, sector, xyz, depe, bmtc);
	
	int n_multi_hits = multi_hit.size()/2;
	
	// closest strip
	//yid[4].id = (int) multi_hit[0];
	yid[3].id = (int) multi_hit[0];
	
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
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
		this_id.id_sharing = multi_hit[3];
		yid.push_back(this_id);
	      }
	    // last id is strip
	    identifier this_id;
	    this_id.name       = yid[3].name;
	    this_id.rule       = yid[3].rule;
	    this_id.id         = (int) multi_hit[2];
	    this_id.time       = yid[3].time;
	    this_id.TimeWindow = yid[3].TimeWindow;
	    this_id.TrackId    = yid[3].TrackId;
	    this_id.id_sharing = multi_hit[3];
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








