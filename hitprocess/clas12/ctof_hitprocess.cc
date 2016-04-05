// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "ctof_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

static ctofConstants initializeCTOFConstants(int runno)
{
	ctofConstants ctc;

	cout<<"Entering initializeCTOF"<<endl;

	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return ctc;

	ctc.runNo      = runno;
	ctc.date       = "2015-11-29";
	if(getenv ("CCDB_CONNECTION") != NULL)
		ctc.connection = (string) getenv("CCDB_CONNECTION");
	else
		ctc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	ctc.variation  = "default";

	ctc.npaddles   = 48;
	ctc.thick      = 3.0;
	
	ctc.dEdxMIP       = 1.956;  // muons in polyvinyltoluene
	ctc.dEMIP         = ctc.thick*ctc.dEdxMIP;
	ctc.pmtPEYld      = 243;
	ctc.pmtQE         = 0.27;
	ctc.pmtDynodeGain = 4.0; 
        ctc.pmtDynodeK    = 0.5; 
	ctc.pmtFactor     = sqrt(1-ctc.pmtQE+(ctc.pmtDynodeK*ctc.pmtDynodeGain+1)/(ctc.pmtDynodeGain-1));
	ctc.tdcLSB        = 41.6667;
	
	cout<<"CTOF:Setting time resolution"<<endl;

	for(int c=1; c<ctc.npaddles+1;c++)
	{
	  ctc.tres.push_back(60.);
	}
	
	
	int isec,ilay,istr;
	
	vector<vector<double> > data;
	
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ctc.connection));
        cout<<"Connecting to "<<ctc.connection<<"/calibration/ctof"<<endl;
	
	cout<<"CTOF:Getting attenuation"<<endl;
	sprintf(ctc.database,"/calibration/ctof/attenuation:%d",ctc.runNo);
	data.clear(); calib->GetCalib(data,ctc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ctc.attlen[isec-1][ilay-1][0].push_back(data[row][3]);
	  ctc.attlen[isec-1][ilay-1][1].push_back(data[row][4]);
	}
	
        cout<<"CTOF:Getting effective_velocity"<<endl;
	sprintf(ctc.database,"/calibration/ctof/effective_velocity:%d",ctc.runNo);
	data.clear(); calib->GetCalib(data,ctc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ctc.veff[isec-1][ilay-1][0].push_back(data[row][3]);
	  ctc.veff[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	cout<<"CTOF:Getting status"<<endl;
	sprintf(ctc.database,"/calibration/ctof/status:%d",ctc.runNo);
	data.clear() ; calib->GetCalib(data,ctc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ctc.status[isec-1][ilay-1][0].push_back(data[row][3]);
	  ctc.status[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	cout<<"CTOF:Getting gain_balance"<<endl;
	sprintf(ctc.database,"/calibration/ctof/gain_balance:%d",ctc.runNo);
	data.clear() ; calib->GetCalib(data,ctc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ctc.countsForMIP[isec-1][ilay-1][0].push_back(data[row][3]);
	  ctc.countsForMIP[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	/* For future use in HitProcess
	cout<<"Getting time_walk"<<endl;
	sprintf(ctc.database,"/calibration/ctof/time_walk:%d",ctc.runNo);
	data.clear() ; calib->GetCalib(data,ctc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ctc.twlk[isec-1][ilay-1][0].push_back(data[row][3]);
	  ctc.twlk[isec-1][ilay-1][1].push_back(data[row][4]);
	  ctc.twlk[isec-1][ilay-1][2].push_back(data[row][5]);
	  ctc.twlk[isec-1][ilay-1][3].push_back(data[row][6]);
	  ctc.twlk[isec-1][ilay-1][4].push_back(data[row][7]);
	  ctc.twlk[isec-1][ilay-1][5].push_back(data[row][8]);
	}
	*/
	
	return ctc;
}

void ctof_HitProcess::initWithRunNumber(int runno)
{
	if(ctc.runNo != runno)
	{
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ctc = initializeCTOFConstants(runno);
		ctc.runNo = runno;
	}
}


map<string, double> ctof_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector = 1;
	int panel  = 1;
	int paddle = identity[0].id;
	
	// Get the paddle length: in ctof paddles are boxes, the length is the y dimension
	double length = aHit->GetDetector().dimensions[2];
	
	trueInfos tInfos(aHit);
	
	// Distances from upstream, downstream
	double dUp = length + tInfos.ly;
	double dDn = length - tInfos.ly;
	
	// attenuation length
	double attlenUp = ctc.attlen[sector-1][panel-1][0][paddle-1];
	double attlenDn = ctc.attlen[sector-1][panel-1][1][paddle-1];

	// attenuation factor
	double attUp  = exp(-dUp/cm/attlenUp);
	double attDn  = exp(-dDn/cm/attlenDn);

	// Attenuated light at PMT
	double eneUp = tInfos.eTot*attUp;
	double eneDn = tInfos.eTot*attDn;

	double adcu = 0.;
	double adcd = 0.;
        double tdcu = 0.;
	double tdcd = 0.;

	// Fluctuate the light measured by the PMT with
	// Poisson distribution for emitted photoelectrons
	// Treat Up and Dn separately, in case nphe=0

	double npheUp = G4Poisson(eneUp*ctc.pmtPEYld);
	eneUp = npheUp/ctc.pmtPEYld;

	if (eneUp>0) {
	                 adcu = eneUp*ctc.countsForMIP[sector-1][panel-1][0][paddle-1]/ctc.dEMIP;
	 //double            A = ctc.twlk[sector-1][panel-1][0][paddle-1];
         //double            B = ctc.twlk[sector-1][panel-1][1][paddle-1];
         //double            C = ctc.twlk[sector-1][panel-1][2][paddle-1];
	 //double   timeWalkUp = A/(B+C*sqrt(adcu));
	  double    timeWalkUp = 0.;
	  double          tUpU = tInfos.time + dUp/ctc.veff[sector-1][panel-1][0][paddle-1]/cm + timeWalkUp;
	  double           tUp = G4RandGauss::shoot(tUpU,  sqrt(2)*ctc.tres[paddle-1]*1e-3);
	                  tdcu = tUp*ctc.tdcLSB;
	}  
	
	double npheDn = G4Poisson(eneDn*ctc.pmtPEYld);
	eneDn = npheDn/ctc.pmtPEYld;

	if (eneDn>0) {
	                 adcd = eneDn*ctc.countsForMIP[sector-1][panel-1][1][paddle-1]/ctc.dEMIP;
	 //double            A = ctc.twlk[sector-1][panel-1][3][paddle-1];
         //double            B = ctc.twlk[sector-1][panel-1][4][paddle-1];
         //double            C = ctc.twlk[sector-1][panel-1][5][paddle-1];
	 //double   timeWalkDn = A/(B+C*sqrt(adcd));
	  double    timeWalkDn = 0.;
	  double          tDnU = tInfos.time + dDn/ctc.veff[sector-1][panel-1][1][paddle-1]/cm + timeWalkDn;
	  double           tDn = G4RandGauss::shoot(tDnU,  sqrt(2)*ctc.tres[paddle-1]*1e-3);
	                  tdcd = tDn*ctc.tdcLSB;
	}  
	
	dgtz["hitn"]   = hitn;
	dgtz["paddle"] = paddle;
	dgtz["ADCL"]   = (int) adcu;
	dgtz["ADCR"]   = (int) adcd;
	dgtz["TDCL"]   = (int) tdcu;
	dgtz["TDCR"]   = (int) tdcd;
	
	return dgtz;
}

vector<identifier>  ctof_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


map< string, vector <int> >  ctof_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// this static function will be loaded first thing by the executable
ctofConstants ctof_HitProcess::ctc = initializeCTOFConstants(-1);
