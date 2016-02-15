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

	int isec,ilay,istr;
	
	// database
	ftc.runNo      = runno;
	ftc.date       = "2015-11-29";
	ftc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	ftc.variation  = "default";

	ftc.npaddles[0] = 23;
	ftc.npaddles[1] = 62;
	ftc.npaddles[2] = 5;

	// Paddle thickness (cm) for MIP peak energy
	ftc.thick[0] = 5.0;
	ftc.thick[1] = 6.0;
	ftc.thick[2] = 5.0; 
	
	vector<vector<double> > data;
	
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ftc.connection));
	
	sprintf(ftc.database,"/calibration/ftof/attenuation:%d",ftc.runNo);
	data.clear(); calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.attlen[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.attlen[isec-1][ilay-1][1].push_back(data[row][4]);
	}
	
	sprintf(ftc.database,"/calibration/ftof/effective_velocity:%d",ftc.runNo);
	data.clear(); calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.veff[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.veff[isec-1][ilay-1][1].push_back(data[row][4]);
	}
	
	sprintf(ftc.database,"/calibration/ftof/status:%d",ftc.runNo);
	data.clear() ; calib->GetCalib(data,ftc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.status[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.status[isec-1][ilay-1][1].push_back(data[row][4]);
	}

	sprintf(ftc.database,"/calibration/ftof/gain_balance:%d",ftc.runNo);
	data.clear() ; calib->GetCalib(data,ftc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  ftc.countsForMIP[isec-1][ilay-1][0].push_back(data[row][3]);
	  ftc.countsForMIP[isec-1][ilay-1][1].push_back(data[row][4]);
	}

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
	

	ftc.dEdxMIP = 1.956 ;  // MeV gm-1 cm-3 (muons in polyvinyltoluene)
	ftc.dEMIP[0] = ftc.thick[0]*ftc.dEdxMIP;
	ftc.dEMIP[1] = ftc.thick[1]*ftc.dEdxMIP;
	ftc.dEMIP[2] = ftc.thick[2]*ftc.dEdxMIP;
	
	// time resolution
	for(int p=0; p<3; p++)
	{
		ftc.sigma0[p] = 0.5 ;
		ftc.sigma1[p] = 0.5 ;
	}
	
	// number of photoelectrons reaching PMT is ~10K / MeV / 4
	ftc.nphePerMevReachingPMT = 2500;
	
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

	
	// Get the paddle length: in ftof paddles are boxes, the length the x
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

	// gain factor to simulate PMT gain matching algorithm
	// i.e.- L,R PMTs are not gain matched to each other, but adjusted so geometeric mean sqrt(L*R)
	// is independent of counter length
	double gainLeft  = sqrt(attLeft*attRight);
	double gainRight = gainLeft;
	
	// multiple of MIP energy attenuated
	double eneL = (tInfos.eTot/ftc.dEMIP[panel-1])*attLeft;
	double eneR = (tInfos.eTot/ftc.dEMIP[panel-1])*attRight;
	
	// Attenuated energy converted to ADC counts
	double adcl = ftc.countsForMIP[sector-1][panel-1][0][paddle-1]*eneL/gainLeft;
	double adcr = ftc.countsForMIP[sector-1][panel-1][1][paddle-1]*eneR/gainRight;

	// timewalk depends on adc values
	double timeWalkLeft  = ftc.twlk[sector-1][panel-1][0][paddle-1]/(1 + ftc.twlk[sector-1][panel-1][1][paddle-1]*sqrt(adcl));
	double timeWalkRight = ftc.twlk[sector-1][panel-1][0][paddle-1]/(1 + ftc.twlk[sector-1][panel-1][1][paddle-1]*sqrt(adcr));
	
	// Unsmeared time at each end using effective velocity and with timewalk correction
	double tLeftU  = tInfos.time/ns + dLeft  / ftc.veff[sector-1][panel-1][0][paddle-1] + timeWalkLeft;
	double tRightU = tInfos.time/ns + dRight / ftc.veff[sector-1][panel-1][1][paddle-1] + timeWalkRight;

	// Fluctuate the light measured by the PMT with Poisson distribution for emitted photoelectrons
	double npheL = G4Poisson(eneL*ftc.nphePerMevReachingPMT);
	double npheR = G4Poisson(eneL*ftc.nphePerMevReachingPMT);
	
	// Smear the time with a gaussian using sigma = A0*A0 + A1*A1/nphe
	double tLeft  = G4RandGauss::shoot(tLeftU,  sqrt(ftc.sigma0[panel-1]*ftc.sigma0[panel-1] + ftc.sigma1[panel-1]*ftc.sigma1[panel-1]/npheL));
	double tRight = G4RandGauss::shoot(tRightU, sqrt(ftc.sigma0[panel-1]*ftc.sigma0[panel-1] + ftc.sigma1[panel-1]*ftc.sigma1[panel-1]/npheR));

	
	double tdclu = tLeftU;
	double tdcru = tRightU;
	double tdcl  = tLeft;
	double tdcr  = tRight;

	
	// applying status
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
	dgtz["ADCLu"]  = adcl;
	dgtz["ADCRu"]  = adcr;
	dgtz["TDCLu"]  = tdclu;
	dgtz["TDCRu"]  = tdcru;
	
	return dgtz;
}

vector<identifier>  ftof_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}




map< string, vector <int> >  ftof_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// this static function will be loaded first thing by the executable
ftofConstants ftof_HitProcess::ftc = initializeFTOFConstants(-1);




