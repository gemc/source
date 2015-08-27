// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"
#include <math.h>

// gemc headers
#include "pcal_hitprocess.h"

static pcConstants initializePCConstants(int runno)
{
	pcConstants pcc;
	pcc.runNo = 0;
	
	pcc.attlen              = 3760.;  // Attenuation Length (mm)
	pcc.TDC_time_to_evio    = 1000.;  // Currently EVIO banks receive time from rol2.c in ps (raw counts x 24 ps/chan. for both V1190/1290), so convert ns to ps.
	pcc.ADC_MeV_to_evio     = 10.  ;  // MIP based calibration is nominally 10 channels/MeV
	pcc.PE_yld              = 11.5 ;  // Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
	pcc.veff                = 160. ;  // Effective velocity of scintillator light (mm/ns)
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
	
	// get sector, stack (inner or outer), view (U, V, W), and strip.
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
	
	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		if(pcc.attlen>0)
		{
			double xlocal = Lpos[s].x();
			if(view==1) latt=pDx2+xlocal;
			if(view==2) latt=pDx2+xlocal;
			if(view==3) latt=pDx2+xlocal;
			//cout<<"xlocal="<<xlocal<<" ylocal="<<ylocal<<" view="<<view<<" strip="<<strip<<" latt="<<latt<<endl;
			Etota = Etota + Edep[s]*exp(-latt/pcc.attlen);
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
	if (tInfos.eTot > 0)
	{
		double PC_npe = G4Poisson(Etota*pcc.PE_yld); //number of photoelectrons
		//  Fluctuations in PMT gain distributed using Gaussian with
		//  sigma SNR = sqrt(ngamma)/sqrt(del/del-1) del = dynode gain = 3 (From RCA PMT Handbook) p. 169)
		//  algorithm, values, and comment above taken from gsim.
		double sigma = sqrt(PC_npe)/1.22;
		double PC_MeV = G4RandGauss::shoot(PC_npe,sigma)*pcc.ADC_MeV_to_evio/pcc.PE_yld;
		if (PC_MeV <= 0) PC_MeV=0.0; // guard against weird, rare events.
		ADC = (int) PC_MeV;
	}
	
	// simulate the tdc.
	TDC = (int) ((tInfos.time+Ttota/tInfos.nsteps)*pcc.TDC_time_to_evio);
	
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
pcConstants pcal_HitProcess::pcc = initializePCConstants(1);







