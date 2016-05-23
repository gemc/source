// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "ft_cal_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> ft_cal_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);

	
	// R.De Vita (April 2009)
	
	// relevant parameter for digitization (in the future should be read from database)
	double tdc_time_to_channel=20;     // conversion factor from time(ns) to TDC channels)
	double tdc_max=8191;               // TDC range
	double time_res=0.2;               // time resolution
	double adc_charge_tochannel=20;    // conversion factor from charge(pC) to ADC channels
	double PbWO4_light_yield =240/MeV; // Lead Tungsten Light Yield (APD have similar QE for
	// fast component, lambda=420nm-ly=120ph/MeV, and slow component,
	// lambda=560nm-ly=20ph/MeV, taking fast component only)
	// double PbWO4_light_yield =672/MeV; // LY at -25 deg=2.8 x LY at +18 deg
	double APD_qe    = 0.70;           // APD Quantum Efficiency (Hamamatsu S8664-55)
	double APD_size  = 100*mm*mm;       // APD size ( 10 mm x 10 mm)
	double APD_gain  = 150;            // based on FT note
	double APD_noise = 0.0033;          // relative noise based on a Voltage and Temperature stability of 10 mV (3.9%/V) and 0.1 C (3.3%/C)
	double AMP_input_noise = 9000;     // preamplifier input noise in number of electrons
	double AMP_gain        = 1800;     // preamplifier gain = 5V/pC x 25 ns (tipical signal duration)
	double light_speed =15;
	
	// Get the crystal length: in the FT crystal are BOXes and the half-length is the 3rd element
	double length = 2 * aHit->GetDetector().dimensions[2];
	// Get the crystal width (rear face): in the FT crystal are BOXes and the half-length is the 2th element
	double width  = 2 * aHit->GetDetector().dimensions[1];
	
	// use Crystal ID to define IDX and IDY
	int IDX = identity[0].id;
	int IDY = identity[1].id;
	
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = 8191;
	
	
	
	if(tInfos.eTot>0)
	{
		/*
		 // commented out to use average time instead of minimum time in TDC calculation
		 for(int s=0; s<nsteps; s++)
		 {
     double dRight = length/2 - Lpos[s].z();              // distance along z between the hit position and the end of the crystal
     double timeR  = times[s] + dRight/cm/light_speed;    // arrival time of the signal at the end of the crystal (speed of light in the crystal=15 cm/ns)
     if(Edep[s]>1*MeV) Tmin=min(Tmin,timeR);              // signal time is set to first hit time with energy above 1 MeV
		 }
		 TDC=int(Tmin*tdc_time_to_channel);
		 if(TDC>tdc_max) TDC=(int)tdc_max;
		 */
		double dRight = length/2 - tInfos.lz;                 // distance along z between the hit position and the end of the crystal
		double timeR  = tInfos.time + dRight/cm/light_speed;  // arrival time of the signal at the end of the crystal (speed of light in the crystal=15 cm/ns)
		// adding spread on time
		timeR=timeR+G4RandGauss::shoot(0., time_res);
		
		TDC=int(timeR*tdc_time_to_channel);
		if(TDC>tdc_max) TDC=(int)tdc_max;
		
		// calculate number of photoelectrons detected by the APD considering the light yield, the q.e., and the size of the sensor
		double npe=G4Poisson(tInfos.eTot*PbWO4_light_yield*0.5*APD_qe*APD_size/width/width);
		//   double npe=(Etot*PbWO4_light_yield*0.5*APD_qe*APD_size/width/width);
		
		// for PMT, an addition factor of 0.5625 is needed to reproduce the 13.5 photoelectrons with a 20% QE
		//   double npe=G4Poisson(Etot*PbWO4_light_yield*0.5*0.5625*APD_qe*APD_size/width/width);
		
		// calculating APD output charge (in number of electrons) and adding noise
		double nel=npe*APD_gain;
		nel=nel*G4RandGauss::shoot(1.,APD_noise);
		if(nel<0) nel=0;
		// adding preamplifier input noise
		nel=nel+AMP_input_noise*G4RandGauss::shoot(0.,1.);
		if(nel<0) nel=0;
		// converting to charge (in picoCoulomb)
		double crg=nel*AMP_gain*1.6e-7;
		// converting to ADC channels
		ADC= (int) (crg*adc_charge_tochannel);
   	
	}
	
	dgtz["hitn"] = hitn;
	dgtz["idx"]  = IDX;
	dgtz["idy"]  = IDY;
	dgtz["adc"]  = ADC;
	dgtz["tdc"]  = TDC;
	
	
	return dgtz;
}

vector<identifier>  ft_cal_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ft_cal_HitProcess :: electronicNoise()
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



map< string, vector <int> >  ft_cal_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}













