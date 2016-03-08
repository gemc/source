// gemc headers
#include "BMT_hitprocess.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"


static bmtConstants initializeBMTConstants(int runno)
{
	// all these constants should be read from CCDB
	bmtConstants bmtc;
	bmtc.runNo = runno;
	
	// all dimensions are in mm
	bmtc.SigmaDrift = 0.4;
	bmtc.hDrift     = 3.0;
	bmtc.hStrip2Det = bmtc.hDrift/2.;
	bmtc.changeFieldScale(-1);  // this needs to be read from DB

	if(runno == -1)
	{
		cout << " > bmt pre-initizialization. " << endl;
		return bmtc;
	}
	
	// Z Detectors
	bmtc.CRZRADIUS[2]	 = 205.8;
	bmtc.CRZNSTRIPS[2] = 768;
	bmtc.CRZSPACING[2] = 0.2;
	bmtc.CRZWIDTH[2]	 = 0.328;
	bmtc.CRZLENGTH[2]  = 444.88;
	bmtc.CRZZMIN[2]	 = -421.75;
	bmtc.CRZZMAX[2]	 = 290.25;
	bmtc.CRZOFFSET[2]	 = 252.1;
	bmtc.CRZXPOS[2]	 = 10.547;
	
	double ZEdge1[bmtc.NREGIONS] = { 30.56*degree, 270.56*degree, 150.56*degree};
	double ZEdge2[bmtc.NREGIONS] = {149.44*degree,  29.44*degree, 269.44*degree};
	
	// Assume the edge boundaries are the same for all regions until get final geometry
	for (int i = 0; i <bmtc.NREGIONS ; ++i)
	{
		for (int j = 0; j <bmtc.NREGIONS ; ++j)
		{
			bmtc.CRZEDGE1[i][j] = ZEdge1[j];
			bmtc.CRZEDGE2[i][j] = ZEdge2[j];
		}
	}
	
	
	// C Detectors
	bmtc.CRCRADIUS[2]	 =	220.8;
	bmtc.CRCNSTRIPS[2] =	1152;
	bmtc.CRCLENGTH[2]	 =	438.6;
	bmtc.CRCSPACING[2] =	0.16;
	bmtc.CRCZMIN[2]	 =	-421.75;
	bmtc.CRCZMAX[2]	 =	290.25;
	bmtc.CRCOFFSET[2]	 =	252.18;

	double CEdge1[bmtc.NREGIONS] = {  30.52*degree, 270.52*degree, 150.52*degree};
	double CEdge2[bmtc.NREGIONS] = { 149.48*degree,  29.48*degree, 269.48*degree};
	
	// Assume the edge boundaries are the same for all regions until get final geometry
	for (int i = 0; i <bmtc.NREGIONS ; ++i)
	{
		for (int j = 0; j <bmtc.NREGIONS ; ++j)
		{
			bmtc.CRCEDGE1[i][j] = CEdge1[j];
			bmtc.CRCEDGE2[i][j] = CEdge2[j];
		}
	}
	
	// pitch CRC --> for the C-detectors, the strips are in bunches of equal pitches
	double CR4C_width[13]={0.345,0.28,0.225,0.175,0.17,0.21,0.26,0.31,0.37,0.44,0.515,0.605,0.7};     // width of the corresponding group of strips
	double CR4C_group[13]={32,32,32,32,624,32,32,32,32,32,32,32,896};                              	  // the number of strips with equal pitches
																																	  // For CR5, no existing value. picked a random value compatible with the geometry
	double CR5C_width[1]={0.253};                                                                     // the number of strips with equal pitches & width of the corresponding group of strips
	double CR5C_group[1]={1024};
	// For CR6 the numbers are final and should not be changed
	double CR6C_width[14]={0.38,0.32,0.27,0.23,0.17,0.18,0.22,0.25,0.29,0.33,0.37,0.41,0.46,0.51};    // the number of strips with equal pitches & width of the corresponding group of strips
	double CR6C_group[14]={32,32,32,32,704,64,32,32,32,32,32,32,32,32};

	int MxGrpSize = 14; // the max number of entries in CRC_group array
	bmtc.CRCGROUP.resize(bmtc.NREGIONS); //3 regions
	bmtc.CRCWIDTH.resize(bmtc.NREGIONS);
	
	for (int i = 0; i <3 ; ++i)
	{
		bmtc.CRCGROUP[i].resize(MxGrpSize);
		bmtc.CRCWIDTH[i].resize(MxGrpSize);
	}
	
	for(int j =0; j<13; j++)
	{
		// region index  0 is CR4
		bmtc.CRCGROUP[0][j] = CR4C_group[j];
		bmtc.CRCWIDTH[0][j] = CR4C_width[j];
	}
	for(int j =0; j<1; j++)
	{
		// region index  1 is CR5
		bmtc.CRCGROUP[1][j] = CR5C_group[j];
		bmtc.CRCWIDTH[1][j] = CR5C_width[j];
	}
	for(int j =0; j<14; j++)
	{
		// region index  2 is CR6
		bmtc.CRCGROUP[2][j] = CR6C_group[j];
		bmtc.CRCWIDTH[2][j] = CR6C_width[j];
	}
	
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
			if(bmtc.runNo == -1)
				cout << " > BMT: Field found with value " << fieldValue[2]/gauss << " gauss. Setting Lorentz angle accordingly." << endl;
			bmtc.changeFieldScale((fieldValue[2]/gauss)/50000.0);
		}
		else
		{
			bmtc.changeFieldScale(0);
			if(bmtc.runNo == -1)
				cout << " > BMT: No field found. Lorentz angle set to zero." << endl;
		}
	}
	
	vector<identifier> yid = id;
	class bmt_strip bmts;
	
	//int layer  = 2*yid[0].id + yid[1].id - 2 ;
	int layer  = yid[0].id;
	int sector = yid[2].id;
	
	double depe = aStep->GetTotalEnergyDeposit();
	//cout << "resolMM " << layer << " " << x << " " << y << " " << z << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;
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




// this static function will be loaded first thing by the executable
bmtConstants BMT_HitProcess::bmtc = initializeBMTConstants(-1);








