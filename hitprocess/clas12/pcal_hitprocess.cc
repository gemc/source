// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"
#include <math.h>

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "pcal_hitprocess.h"

static pcConstants initializePCConstants(int runno)
{
	pcConstants pcc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return pcc;
		
	int isec,ilay,istr;
	
	// database
	pcc.runNo      = runno;
	pcc.date       = "2015-11-29";
	if(getenv ("CCDB_CONNECTION") != NULL)
		pcc.connection = (string) getenv("CCDB_CONNECTION");
	else
		pcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	pcc.variation  = "default";
	
	pcc.TDC_time_to_evio    = 1000.;  // Currently EVIO banks receive time from rol2.c in ps (raw counts x 24 ps/chan. for both V1190/1290), so convert ns to ps.
	pcc.ADC_MeV_to_evio     = 10.  ;  // MIP based calibration is nominally 10 channels/MeV
	pcc.veff                = 160. ;  // Effective velocity of scintillator light (mm/ns)
	pcc.pmtPEYld            = 11.5 ; // Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
	pcc.pmtQE               = 0.27 ; 
	pcc.pmtDynodeGain       = 4.0  ; 
        pcc.pmtDynodeK          = 0.5  ; // K=0 (Poisson) K=1(exponential)	vector<vector<double> > data;
	//  Fluctuations in PMT gain distributed using Gaussian with
	//  sigma=1/SNR where SNR = sqrt[(1-QE+(k*del+1)/(del-1))/npe] del = dynode gain k=0-1
	//  Adapted from G-75 (pg. 169) and and G-111 (pg. 174) from RCA PMT Handbook.
	//  Factor k for dynode statistics can range from k=0 (Poisson) to k=1 (exponential).
	//  Note: GSIM sigma was incorrect (used 1/sigma for sigma).
	pcc.pmtFactor           = sqrt(1-pcc.pmtQE+(pcc.pmtDynodeK*pcc.pmtDynodeGain+1)/(pcc.pmtDynodeGain-1));

	vector<vector<double> > data;
        auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(pcc.connection));
	
	sprintf(pcc.database,"/calibration/ec/attenuation:%d",pcc.runNo);
	data.clear(); calib->GetCalib(data,pcc.database);
	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  pcc.attlen[isec-1][ilay-1][0].push_back(data[row][3]);
	  pcc.attlen[isec-1][ilay-1][1].push_back(data[row][4]);
	  pcc.attlen[isec-1][ilay-1][2].push_back(data[row][5]);
        }
	
	return pcc;
}


void pcal_HitProcess::initWithRunNumber(int runno)
{
	if(pcc.runNo != runno)
	{
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		pcc = initializePCConstants(runno);
		pcc.runNo = runno;
	}
}
map<string, double> pcal_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// get sector, view (U, V, W), and strip.
	int sector = identity[0].id;
	int module = identity[1].id;
	int view   = identity[2].id;
	int strip  = identity[3].id;
	trueInfos tInfos(aHit);
	
	HCname = "PCAL Hit Process";
	
	// Get scintillator volume x dimension (mm)
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.
	
	//cout<<"pDx2="<<pDx2<<" pDy1="<<pDy1<<endl;
	
	
	// Get Total Energy deposited
	double Etota = 0;
	double Ttota = 0;
	double latt = 0;
	
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	
	double att;
	
	double A = pcc.attlen[sector-1][view-1][0][strip-1];
	double B = pcc.attlen[sector-1][view-1][1][strip-1]*10.;
	double C = pcc.attlen[sector-1][view-1][2][strip-1];

	//cout<<"sector "<<sector<<"view "<<view<<"strip "<<strip<<"B "<<B<<endl;
	
	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		if(B>0)
		{
			double xlocal = Lpos[s].x();
			if(view==1) latt=pDx2+xlocal;
			if(view==2) latt=pDx2+xlocal;
			if(view==3) latt=pDx2+xlocal;
			att   = A*exp(-latt/B)+C;
			Etota = Etota + Edep[s]*att;
			Ttota = Ttota + latt/pcc.veff;
		}
		else
		{
			Etota = Etota + Edep[s];
		}
	}
	
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = 0;
	
	// simulate the adc value.
	if (Etota > 0) {
	  double PC_npe = G4Poisson(Etota*pcc.pmtPEYld); //number of photoelectrons
	  if (PC_npe>0) {
	    double sigma = pcc.pmtFactor/sqrt(PC_npe);
	    double PC_MeV = G4RandGauss::shoot(PC_npe,sigma)*pcc.ADC_MeV_to_evio/pcc.pmtPEYld;
	    if (PC_MeV>0) {
	      ADC = (int) PC_MeV;
	      TDC = (int) ((tInfos.time+Ttota/tInfos.nsteps)*pcc.TDC_time_to_evio);
	    }
	  }
	}
	
	// EVIO banks record time with offset determined by position of data in capture window.  On forward carriage this is currently
	// around 1.4 us.  This offset is omitted in the simulation.
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["module"] = module;
	dgtz["view"]   = view;
	dgtz["strip"]  = strip;
	dgtz["ADC"]    = ADC;
	dgtz["TDC"]    = TDC;
	
	//cout << "sector = " << sector << " layer = " << module << " view = " << view << " strip = " << strip << " PL_ADC = " << ADC << " TDC = " << TDC << " Edep = " << Etot << endl;
	
	return dgtz;
}


vector<identifier>  pcal_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	//int sector             = yid[0].id;
	//int layer              = yid[1].id;
	//int view               = yid[2].id; // get the view: U->1, V->2, W->3
	//int strip	       = yid[3].id;
	//return yid;
	id[id.size()-1].id_sharing = 1;
	return id;
}



map< string, vector <int> >  pcal_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// this static function will be loaded first thing by the executable
pcConstants pcal_HitProcess::pcc = initializePCConstants(-1);







