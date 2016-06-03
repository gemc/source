// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "ftof_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

static ftofConstants initializeFTOFConstants(int runno)
{
	ftofConstants ftc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return ftc;

	ftc.runNo      = runno;
	ftc.date       = "2015-11-29";
	if(getenv ("CCDB_CONNECTION") != NULL)
		ftc.connection = (string) getenv("CCDB_CONNECTION");
	else
		ftc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	ftc.variation  = "default";

	ftc.npaddles[0] = 23;
	ftc.npaddles[1] = 62;
	ftc.npaddles[2] = 5;

	ftc.thick[0] = 5.0;
	ftc.thick[1] = 6.0;
	ftc.thick[2] = 5.0; 
	
	ftc.dEdxMIP       = 1.956;  // muons in polyvinyltoluene
	ftc.pmtPEYld      = 500;
	ftc.tdcLSB        = 41.6667; // counts per ns (24 ps LSB)
	
	cout<<"FTOF:Setting time resolution"<<endl;
	for(int p=0; p<3; p++)
	{
	  for(int c=1; c<ftc.npaddles[p]+1;c++)
	    {
	      if(p==0) ftc.tres[p].push_back(1e-3*(c*5.45+74.55)); //ps to ns
	      if(p==1) ftc.tres[p].push_back(1e-3*(c*0.90+29.10)); //ps to ns
	      if(p==2) ftc.tres[p].push_back(1e-3*(c*5.00+145.0)); //ps to ns
	    }
	}
	
	ftc.dEMIP[0] = ftc.thick[0]*ftc.dEdxMIP;
	ftc.dEMIP[1] = ftc.thick[1]*ftc.dEdxMIP;
	ftc.dEMIP[2] = ftc.thick[2]*ftc.dEdxMIP;
	
	int isec,ilay,istr;
	
	vector<vector<double> > data;
	
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ftc.connection));
        cout<<"Connecting to "<<ftc.connection<<"/calibration/ftof"<<endl;

	cout<<"FTOF:Getting attenuation"<<endl;
	sprintf(ftc.database,"/calibration/ftof/attenuation:%d",ftc.runNo);
	data.clear(); calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.attlen[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.attlen[isec-1][ilay-1][1].push_back(data[row][4]);
	}
	
        cout<<"FTOF:Getting effective_velocity"<<endl;
	sprintf(ftc.database,"/calibration/ftof/effective_velocity:%d",ftc.runNo);
	data.clear(); calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.veff[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.veff[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	cout<<"FTOF:Getting status"<<endl;
	sprintf(ftc.database,"/calibration/ftof/status:%d",ftc.runNo);
	data.clear() ; calib->GetCalib(data,ftc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.status[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.status[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	cout<<"FTOF:Getting gain_balance"<<endl;
	sprintf(ftc.database,"/calibration/ftof/gain_balance:%d",ftc.runNo);
	data.clear() ; calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.countsForMIP[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.countsForMIP[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	cout<<"FTOF:Getting time_walk"<<endl;
	sprintf(ftc.database,"/calibration/ftof/time_walk:%d",ftc.runNo);
	data.clear() ; calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.twlk[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.twlk[isec-1][ilay-1][1].push_back(data[row][4]);
	  ftc.twlk[isec-1][ilay-1][2].push_back(data[row][5]);
	  ftc.twlk[isec-1][ilay-1][3].push_back(data[row][6]);
	  ftc.twlk[isec-1][ilay-1][4].push_back(data[row][7]);
	  ftc.twlk[isec-1][ilay-1][5].push_back(data[row][8]);
	}
	
	return ftc;
}

void ftof_HitProcess::initWithRunNumber(int runno)
{
	if(ftc.runNo != runno)
	{
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ftc = initializeFTOFConstants(runno);
		ftc.runNo = runno;
	}
}


map<string, double> ftof_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{	
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	int sector = identity[0].id;
	int panel  = identity[1].id;
	int paddle = identity[2].id;
	trueInfos tInfos(aHit);
	
	// Get the paddle half-length 
	double length = aHit->GetDetector().dimensions[0];
	
	// Distances from left, right
	double dLeft  = length + tInfos.lx;
	double dRight = length - tInfos.lx;
		
	// attenuation length
	double attlenL = ftc.attlen[sector-1][panel-1][0][paddle-1];
	double attlenR = ftc.attlen[sector-1][panel-1][1][paddle-1];
	
	// attenuation factor
	double attLeft  = exp(-dLeft/cm/attlenL);
	double attRight = exp(-dRight/cm/attlenR);

	// Gain factors to simulate FTOF PMT gain matching algorithm.
	// Each L,R PMT pair has HV adjusted so geometeric mean sqrt(L*R)
	// is independent of counter length, which compensates for
	// the factor exp(-L/2/attlen) where L=full length of bar.
	double gainLeft  = sqrt(attLeft*attRight);
	double gainRight = gainLeft;

	// Attenuated light at PMT
	double eneL = tInfos.eTot*attLeft;
	double eneR = tInfos.eTot*attRight;

	// giving geantinos some energies
	if(aHit->GetPID() == 0)
	{
		double gmomentum = aHit->GetMom().mag()/GeV;
		eneL = gmomentum*attLeft;
		eneR = gmomentum*attRight;

	}
		
	double adcl  = 0;
	double adcr  = 0;
	double adclu = 0;
	double adcru = 0;
	double tdcl  = 0;
	double tdcr  = 0;
	double tdclu = 0;
	double tdcru = 0;
	
	// Fluctuate the light measured by the PMT with
	// Poisson distribution for emitted photoelectrons
	// Treat L and R separately, in case nphe=0

	if (eneL>0)
		adclu = eneL*ftc.countsForMIP[sector-1][panel-1][0][paddle-1]/ftc.dEMIP[panel-1]/gainLeft;
	
	if (eneR>0)
		adcru = eneR*ftc.countsForMIP[sector-1][panel-1][1][paddle-1]/ftc.dEMIP[panel-1]/gainRight;
	
	
	double npheL = G4Poisson(eneL*ftc.pmtPEYld);
	eneL = npheL/ftc.pmtPEYld;

        if (eneL>0) {
	                 adcl = eneL*ftc.countsForMIP[sector-1][panel-1][0][paddle-1]/ftc.dEMIP[panel-1]/gainLeft;
	  double            A = ftc.twlk[sector-1][panel-1][0][paddle-1];
	  double            B = ftc.twlk[sector-1][panel-1][1][paddle-1];
	  //double            C = ftc.twlk[sector-1][panel-1][2][paddle-1];
	  double timeWalkLeft = A/pow(adcl,B);
	  double       tLeftU = tInfos.time + dLeft/ftc.veff[sector-1][panel-1][0][paddle-1]/cm + timeWalkLeft;
	  double        tLeft = G4RandGauss::shoot(tLeftU,  sqrt(2)*ftc.tres[panel-1][paddle-1]);
	                tdclu = tLeftU*ftc.tdcLSB;
	                 tdcl = tLeft*ftc.tdcLSB;
	}  
	
	double npheR = G4Poisson(eneR*ftc.pmtPEYld);
	eneR = npheR/ftc.pmtPEYld;

	if (eneR>0) {
	                  adcr = eneR*ftc.countsForMIP[sector-1][panel-1][1][paddle-1]/ftc.dEMIP[panel-1]/gainRight;
	  double             A = ftc.twlk[sector-1][panel-1][3][paddle-1];
	  double             B = ftc.twlk[sector-1][panel-1][4][paddle-1];
	  //double             C = ftc.twlk[sector-1][panel-1][5][paddle-1];	
	  double timeWalkRight = A/pow(adcr,B);	
	  double       tRightU = tInfos.time + dRight/ftc.veff[sector-1][panel-1][1][paddle-1]/cm + timeWalkRight;	
	  double        tRight = G4RandGauss::shoot(tRightU, sqrt(2)*ftc.tres[panel-1][paddle-1]);	
	                 tdcru = tRightU*ftc.tdcLSB;
	                  tdcr = tRight*ftc.tdcLSB;
	}
	
	// Status flags
	switch (ftc.status[sector-1][panel-1][0][paddle-1])
	{
		case 0:
			break;
		case 1:
			adcl = 0;
			break;
		case 2:
			tdcl = 0;
			break;
		case 3:
			adcl = tdcl = 0;
			break;

		case 5:
			break;
			
		default:
			cout << " > Unknown FTOF status: " << ftc.status[sector-1][panel-1][0][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " left " << endl;
	}
	
	switch (ftc.status[sector-1][panel-1][1][paddle-1])
	{
		case 0:
			break;
		case 1:
			adcr = 0;
			break;
		case 2:
			tdcr = 0;
			break;
		case 3:
			adcr = tdcr = 0;
			break;
			
		case 5:
			break;
			
		default:
			cout << " > Unknown FTOF status: " << ftc.status[sector-1][panel-1][1][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " right " << endl;
	}
	
//	cout << " > FTOF status: " << ftc.status[sector-1][panel-1][0][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " left: " << adcl << endl;
//	cout << " > FTOF status: " << ftc.status[sector-1][panel-1][1][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " right:  " << adcr << endl;
	
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["paddle"] = paddle;
	dgtz["ADCL"]   = adcl;
	dgtz["ADCR"]   = adcr;
	dgtz["TDCL"]   = tdcl;
	dgtz["TDCR"]   = tdcr;
	dgtz["ADCLu"]  = adclu;
	dgtz["ADCRu"]  = adcru;
	dgtz["TDCLu"]  = tdclu;
	dgtz["TDCRu"]  = tdcru;
	
	return dgtz;
}

vector<identifier>  ftof_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ftof_HitProcess :: electronicNoise()
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



map< string, vector <int> >  ftof_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// this static function will be loaded first thing by the executable
ftofConstants ftof_HitProcess::ftc = initializeFTOFConstants(-1);




