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
	
	ctc.dEdxMIP       = 1.956;   // muons in polyvinyltoluene
	ctc.dEMIP         = ctc.thick*ctc.dEdxMIP;
	ctc.pmtPEYld      = 500;
	ctc.tdcLSB        = 41.6667; // counts per ns (24 ps LSB)
	
	cout<<"CTOF:Setting time resolution"<<endl;

	for(int c=1; c<ctc.npaddles+1;c++)
	{
	  ctc.tres.push_back(1e-3*65.); //ps to ns
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

	ctc.lengthHighPitch = 35.013*25.4/2;  // length of long bar
	ctc.lengthLowPitch  = 34.664*25.4/2;  // length of short bar

	// setting voltage signal parameters
	ctc.vpar[0] = 50;  // delay, ns
	ctc.vpar[1] = 10;  // rise time, ns
	ctc.vpar[2] = 20;  // fall time, ns
	ctc.vpar[3] = 1;   // amplifier

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

	// odd numbered paddles are short
	// even numbered are long
	double length = ctc.lengthLowPitch;
	if(paddle%2 == 0) length = ctc.lengthHighPitch;

//
//	if(aHit->isElectronicNoise)
//	{
//		dgtz["hitn"]   = -hitn;
//		dgtz["paddle"] = paddle;
//
//		dgtz["ADCU"]   = (int) aHit->GetEdep()[0]*100;
////		dgtz["ADCD"]   = (int) aHit->GetE()*100;
////		dgtz["TDCU"]   = (int) aHit->GetE()*100;
////		dgtz["TDCD"]   = (int) aHit->GetE()*100;
////		dgtz["ADCUu"]  = (int) aHit->GetE()*100;
////		dgtz["ADCDu"]  = (int) aHit->GetE()*100;
////		dgtz["TDCUu"]  = (int) aHit->GetE()*100;
////		dgtz["TDCDu"]  = (int) aHit->GetE()*100;
//
//		return dgtz;
//	}

	// Get the paddle length: in ctof paddles are boxes, the length is the y dimension
	// double length = aHit->GetDetector().dimensions[2];




	trueInfos tInfos(aHit);

	// Distances from upstream, downstream
	// ctof paddle center is exactly the target position,
	// so z is also the local coordinate
	double dUp = length + tInfos.z;
	double dDn = length - tInfos.z;
	
	// attenuation length
	double attlenUp = ctc.attlen[sector-1][panel-1][0][paddle-1];
	double attlenDn = ctc.attlen[sector-1][panel-1][1][paddle-1];

	// attenuation factor
	double attUp  = exp(-dUp/cm/attlenUp);
	double attDn  = exp(-dDn/cm/attlenDn);

	// Gain factors to simulate CTOF PMT gain matching algorithm.
	// Each U,D PMT pair has HV adjusted so geometeric mean sqrt(U*D)
	// is independent of counter length, which compensates for
	// the factor exp(-L/2/attlen) where L=full length of bar.
	double gainUp = sqrt(attUp*attDn);
	double gainDn = gainUp;

	// Attenuated light at PMT
	double eneUp = tInfos.eTot*attUp;
	double eneDn = tInfos.eTot*attDn;

	double adcu = 0.;
	double adcd = 0.;
        double tdcu = 0.;
	double tdcd = 0.;
	double adcuu = 0.;
	double adcdu = 0.;
	double tdcuu = 0.;
	double tdcdu = 0.;

	// Fluctuate the light measured by the PMT with
	// Poisson distribution for emitted photoelectrons
	// Treat Up and Dn separately, in case nphe=0

	if (eneUp>0)
		adcuu = eneUp*ctc.countsForMIP[sector-1][panel-1][0][paddle-1]/ctc.dEMIP/gainUp;
	
	if (eneDn>0)
		adcdu = eneDn*ctc.countsForMIP[sector-1][panel-1][1][paddle-1]/ctc.dEMIP/gainDn;


	double npheUp = G4Poisson(eneUp*ctc.pmtPEYld);
	eneUp = npheUp/ctc.pmtPEYld;

	if (eneUp>0) {
	                 adcu = eneUp*ctc.countsForMIP[sector-1][panel-1][0][paddle-1]/ctc.dEMIP/gainUp;
	 //double            A = ctc.twlk[sector-1][panel-1][0][paddle-1];
         //double            B = ctc.twlk[sector-1][panel-1][1][paddle-1];
         //double            C = ctc.twlk[sector-1][panel-1][2][paddle-1];
	 //double   timeWalkUp = A/(B+C*sqrt(adcu));
	  double    timeWalkUp = 0.;
	  double          tUpU = tInfos.time + dUp/ctc.veff[sector-1][panel-1][0][paddle-1]/cm + timeWalkUp;
	  double           tUp = G4RandGauss::shoot(tUpU,  sqrt(2)*ctc.tres[paddle-1]);
	                 tdcuu = tUpU*ctc.tdcLSB;
	                  tdcu = tUp*ctc.tdcLSB;
	}  
	
	double npheDn = G4Poisson(eneDn*ctc.pmtPEYld);
	eneDn = npheDn/ctc.pmtPEYld;

	if (eneDn>0) {
	                 adcd = eneDn*ctc.countsForMIP[sector-1][panel-1][1][paddle-1]/ctc.dEMIP/gainDn;
	 //double            A = ctc.twlk[sector-1][panel-1][3][paddle-1];
         //double            B = ctc.twlk[sector-1][panel-1][4][paddle-1];
         //double            C = ctc.twlk[sector-1][panel-1][5][paddle-1];
	 //double   timeWalkDn = A/(B+C*sqrt(adcd));
	  double    timeWalkDn = 0.;
	  double          tDnU = tInfos.time + dDn/ctc.veff[sector-1][panel-1][1][paddle-1]/cm + timeWalkDn;
	  double           tDn = G4RandGauss::shoot(tDnU,  sqrt(2)*ctc.tres[paddle-1]);
	                 tdcdu = tDnU*ctc.tdcLSB;
	                  tdcd = tDn*ctc.tdcLSB;
	}  
	
	dgtz["hitn"]   = hitn;
	dgtz["paddle"] = paddle;
	dgtz["ADCU"]   = (int) adcu;
	dgtz["ADCD"]   = (int) adcd;
	dgtz["TDCU"]   = (int) tdcu;
	dgtz["TDCD"]   = (int) tdcd;
	dgtz["ADCUu"]  = (int) adcuu;
	dgtz["ADCDu"]  = (int) adcdu;
	dgtz["TDCUu"]  = (int) tdcuu;
	dgtz["TDCDu"]  = (int) tdcdu;
	
	return dgtz;
}

vector<identifier>  ctof_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ctof_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;

	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)

//	for(unsigned int p=1; p<10; p++)
//	{
//		vector<identifier> thisID;
//
//		// for paddle, identifier is only 1 dimensional: paddle ID
//		identifier thisIdentifier;
//		thisIdentifier.id = p;
////		identifier.name = "ctofNoise";
//
//		thisID.push_back(thisIdentifier);
//
//		MHit* thisNoiseHit = new MHit(10.0*p, (double) p, thisID, p);
//
//		noiseHits.push_back(thisNoiseHit);
//	}
//
	return noiseHits;
}

map< string, vector <int> >  ctof_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > ctof_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double ctof_HitProcess :: voltage(double charge, double time, double forTime)
{
	//	return 0.0;
	return DGauss(forTime, ctc.vpar, charge, time);
}

// this static function will be loaded first thing by the executable
ctofConstants ctof_HitProcess::ctc = initializeCTOFConstants(-1);




