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
	
	// database
	ftc.runNo = runno;
	ftc.date       = "2015-11-15";
	ftc.connection = "mysql://clas12writer:geom3try@clasdb.jlab.org/clas12";
	ftc.variation  = "main";

	// temporary reading attenuations, these are dummy numbers and don't set the values yet
	ftc.database   = "/calibration/ftof/attenuation";

	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ftc.connection));

	
	
	
	auto_ptr<Assignment> calibModel(calib->GetAssignment(ftc.database));
	
	
	//table holds information about columns
//	for(size_t rowI = 0; rowI < calibModel->GetRowsCount(); rowI++)
//	{
//		cout << "  sector: "       << calibModel->GetValueInt(rowI, 0)
//		     << "  panel:  "       << calibModel->GetValue(rowI, 1)
//		     << "  paddle:  "      << calibModel->GetValueInt(rowI, 2)
//			  << "  length_left:  " << calibModel->GetValueDouble(rowI, 3) << endl;
//	}
//	

	
	
	ftc.npaddles[0] = 23;
	ftc.npaddles[1] = 62;
	ftc.npaddles[2] = 5;
	
	for(int s=0; s<6; s++)
		for(int p=0; p<3; p++)
			for(int lr=0; lr<2; lr++)
				for(int i=0; i<ftc.npaddles[p]; i++)
					ftc.status[s][p][lr].push_back(0);
	
	
	if(runno == 13)
	{
		// checking status: assigning 1 to sector 1 panel 1 left paddle 11
		ftc.status[0][0][0][10] = 3;
		ftc.status[0][0][1][10] = 3;
	}
	if(runno == 22)
	{
		// checking status: assigning 1 to sector 1 panel 1 left paddle 13
		ftc.status[0][0][0][12] = 3;
		ftc.status[0][0][1][12] = 3;
	}
	if(runno == 30)
	{
		// checking status: assigning 1 to sector 1 panel 1 left paddle 15
		ftc.status[0][0][0][14] = 3;
		ftc.status[0][0][1][14] = 3;
	}
	
	// effective velocity
	for(int s=0; s<6; s++)
		for(int p=0; p<3; p++)
			for(int i=0; i<ftc.npaddles[p]; i++)
				ftc.veff[s][p].push_back(16);
	
	// counts to minimum ionizing
	for(int s=0; s<6; s++)
		for(int p=0; p<3; p++)
			for(int lr=0; lr<2; lr++)
				for(int i=0; i<ftc.npaddles[p]; i++)
					ftc.countsForAMinimumIonizing[s][p][lr].push_back(2000);
	
	// attenuation length
	for(int p=0; p<3; p++)
	{
		ftc.attLengthPars[0][p] = 81.725 ;
		ftc.attLengthPars[1][p] = 0.35771 ;
	}
	
	// de/dx = 2MeV / g/ cm3 for MIP in the FTOF scintillators
	ftc.dEdxMIP = 2;
	
	// time walk correction
	// counts to minimum ionizing
	for(int s=0; s<6; s++)
		for(int p=0; p<3; p++)
			for(int lr=0; lr<2; lr++)
				for(int i=0; i<ftc.npaddles[p]; i++)
				{
					ftc.twlk_A0[s][p][lr].push_back(50.0);
					ftc.twlk_A1[s][p][lr].push_back(0.852);
				}
	
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
	double attLength = ftc.attLengthPars[0][panel-1] + length/cm*ftc.attLengthPars[1][panel-1];
	
	// attenuation factor
	double attLeft  = exp(-dLeft/cm/attLength);
	double attRight = exp(-dRight/cm/attLength);

	// gain factor to simulate PMT gain matching algorithm
	// i.e.- L,R PMTs are not gain matched to each other, but adjusted so geometeric mean sqrt(L*R)
	// is independent of counter length
	double gainLeft  = sqrt(attLeft*attRight);
	double gainRight = gainLeft;
	
	// multiple of MIP energy attenuated
	double eneL = (tInfos.eTot/ftc.dEdxMIP)*attLeft;
	double eneR = (tInfos.eTot/ftc.dEdxMIP)*attRight;
	
	// attenuated energy, converted in counts
	double adcl = ftc.countsForAMinimumIonizing[sector-1][panel-1][0][paddle-1]*eneL/gainLeft;
	double adcr = ftc.countsForAMinimumIonizing[sector-1][panel-1][1][paddle-1]*eneR/gainRight;

	// timewalk depends on adc values
	double timeWalkLeft  = ftc.twlk_A0[sector-1][panel-1][0][paddle-1]/(1 + ftc.twlk_A1[sector-1][panel-1][0][paddle-1]*sqrt(adcl));
	double timeWalkRight = ftc.twlk_A0[sector-1][panel-1][1][paddle-1]/(1 + ftc.twlk_A1[sector-1][panel-1][1][paddle-1]*sqrt(adcr));
	
	// time at each end using effective velocity and with timewalk correction
	// (Unsmeared)
	double tLeftU  = tInfos.time/ns + dLeft  / ftc.veff[sector-1][panel-1][paddle-1] + timeWalkLeft;
	double tRightU = tInfos.time/ns + dRight / ftc.veff[sector-1][panel-1][paddle-1] + timeWalkRight;
	
	double npheL = G4Poisson(eneL*ftc.nphePerMevReachingPMT);
	double npheR = G4Poisson(eneL*ftc.nphePerMevReachingPMT);

	
	// smearing the time by a gaussian with
	// sigma = A0*A0 + A1*A1/nphe
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
ftofConstants ftof_HitProcess::ftc = initializeFTOFConstants(1);




