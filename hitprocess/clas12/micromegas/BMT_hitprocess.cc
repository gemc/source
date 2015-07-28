// gemc headers
#include "BMT_hitprocess.h"
#include "bmt_strip.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

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
	
	return dgtz;
}



vector<identifier>  BMT_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	double x, y, z;
	G4ThreeVector   xyz    = aStep->GetPostStepPoint()->GetPosition();
	x = xyz.x()/mm;
	y = xyz.y()/mm;
	z = xyz.z()/mm;
	
	vector<identifier> yid = id;
	class bmt_strip bmts;
	bmts.fill_infos();
	
	//int layer  = 2*yid[0].id + yid[1].id - 2 ;
	int layer  = yid[0].id;
	int sector = yid[2].id;
	
	//yid[3].id = bmts.FindStrip(layer-1, sector-1, x, y, z);
	double depe = aStep->GetTotalEnergyDeposit();
	//cout << "resolMM " << layer << " " << x << " " << y << " " << z << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;
	vector<double> multi_hit = bmts.FindStrip(layer-1, sector-1, x, y, z, depe);
	
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



map< string, vector <int> >  BMT_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}












