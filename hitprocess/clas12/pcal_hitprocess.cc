// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"
#include <math.h>

// gemc headers
#include "pcal_hitprocess.h"

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
	
	// Attenuation Length (mm)
	double attlen=3760.;
	
	// Get scintillator volume x dimension (mm)
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.
	
	//cout<<"pDx2="<<pDx2<<" pDy1="<<pDy1<<endl;
	
	
	// Get Total Energy deposited
	double Etota = 0;
	double latt = 0;
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	
	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		if(attlen>0)
		{
			double xlocal = Lpos[s].x();
			if(view==1) latt=pDx2+xlocal;
			if(view==2) latt=pDx2+xlocal;
			if(view==3) latt=pDx2+xlocal;
			//cout<<"xlocal="<<xlocal<<" ylocal="<<ylocal<<" view="<<view<<" strip="<<strip<<" latt="<<latt<<endl;
			Etota = Etota + Edep[s]*exp(-latt/attlen);
		}
		else
		{
			Etota = Etota + Edep[s];
		}
	}
	
	// adapted by Alex Piaseczny, Canisius College Medium Energy Nuclear Physics (CMENP), through
	// Dr. Michael Wood from Gilfoyle EC code
	// Jerry Gilfoyle, Feb, 2010
	
	// parameters: factor  - conversion for adc (MeV/channel).
	//             pcal_tdc_to_channel - conversion factor for tdc (ns/channel).
	
	// int pcal_tdc_time_to_channel = (int) gpars["PCAL/pcal_tdc_time_to_channel"];     // conversion from time (ns) to TDC channels. name has to match parameters
	
	double TDC_time_to_channel = 20 ;
	double PCfactor = 11.5;               // number of p.e. divided by the energy deposited in MeV; measured from PCAL cosmic tests
	int TDC_MAX = 4095;                // max value for PC tdc.
	double pc_MeV_to_channel = 10.;       // conversion from energy (MeV) to ADC channels
	
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = TDC_MAX;
	
	// simulate the adc value.
	if (tInfos.eTot > 0)
	{
		double PC_npe = G4Poisson(Etota*PCfactor); //number of photoelectrons
		//  Fluctuations in PMT gain distributed using Gaussian with
		//  sigma SNR = sqrt(ngamma)/sqrt(del/del-1) del = dynode gain = 3 (From RCA PMT Handbook) p. 169)
		//  algorithm, values, and comment above taken from gsim.
		double sigma = sqrt(PC_npe)/1.22;
		double PC_MeV = G4RandGauss::shoot(PC_npe,sigma)*pc_MeV_to_channel/PCfactor;
		if (PC_MeV <= 0) PC_MeV=0.0; // guard against weird, rare events.
		ADC = (int) PC_MeV;
	}
	
	// simulate the tdc.
	TDC = (int) (tInfos.time*TDC_time_to_channel);
	if (TDC > TDC_MAX) TDC = TDC_MAX;
	
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











