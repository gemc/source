// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "cormo_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> cormo_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{ 
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector = identity[0].id;
	int layer  = identity[1].id;
	int paddle = identity[2].id;
		
	// Digitization Parameters
	double light_yield=1600/MeV;             // number of optical photons pruced in the scintillator per MeV of deposited energy
	double att_length=108*cm;                 // light at tenuation length
	double sensor_surface=pow(3.5*cm,2)*pi;   // area of photo sensor
	double paddle_surface=pow(10*cm,2);       // paddle surface
	double light_coll=sensor_surface/paddle_surface;  // ratio of photo_sensor area over paddle section ~ light collection efficiency

	double sensor_qe=0.15;                     // photo sensor quantum efficiency
	double sensor_gain=0.472*1.6*1.275;         // pmt gain x electron charge in pC (4x10^6)x(1.6x10^-7) ->val * 1.6 * 0.6 = 0,963   modifica 21/01/10 
	double adc_conv=10;                       // conversion factor from pC to ADC (typical sensitivy of CAEN VME QDC is of 0.1 pC/ch)
	double adc_ped=0;                         // ADC Pedestal
	double veff=13*cm/ns;                     // light velocity in scintillator
	double tdc_conv=1/0.001/ns;               // TDC conversion factor
	
	
	// initialize ADC and TDC
	double etotL = 0;
	double etotR = 0;
	double timeL = 0;
	double timeR = 0;
	int ADCL = 0;
	int ADCR = 0;
	int TDCL = 4096;
	int TDCR = 4096;
	
	double etotB = 0; // signal propagating to backward end of paddle hit happened in
	double etotF = 0; // signal propagating to forward end of the paddle the hit happened in, round u-turn and along neighbouring paddle
	double timeB = 0;
	double timeF = 0;
	int ADCB = 0;
	int ADCF = 0;
	int TDCB = 4096;
	int TDCF = 4096;
	
	
	// Get the paddle length: in cormo paddles are along y
//	double length = aHit->GetDetector().dimensions[2];
	double length = 20*cm;
	
	// Get info about detector material to eveluate Birks effect
	double birks_constant=aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
	//	cout << " Birks constant is: " << birks_constant << endl;
	//	cout << aHit->GetDetector().GetLogical()->GetMaterial()->GetName() << endl;
	
	
	double time_min[4] = {0,0,0,0};
	
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double>      Edep = aHit->GetEdep();
	vector<G4double>      Dx   = aHit->GetDx();
	// Charge for each step
	vector<int> charge = aHit->GetCharges();
	vector<G4double> times = aHit->GetTime();
	
	unsigned int nsteps = Edep.size();
	double       Etot   = 0;
	for(unsigned int s=0; s<nsteps; s++) Etot = Etot + Edep[s];

	if(Etot>0)
	{
	  for(unsigned int s=0; s<nsteps; s++)
		{
			// Distances from left, right
			double dLeft  = length + Lpos[s].y();
			double dRight = length - Lpos[s].y();
			
			// cout << "\n Distances: " << endl;
			// cout << "\t dLeft, dRight, dBack, dFwd: " << dLeft << ", " << dRight << ", " << dBack << ", " << dFwd << endl;
			
			// apply Birks effect
// 			double stepl = 0.;
			
//			if (s == 0){
//				stepl = sqrt(pow((Lpos[s+1].x() - Lpos[s].x()),2) + pow((Lpos[s+1].y() - Lpos[s].y()),2) + pow((Lpos[s+1].z() - Lpos[s].z()),2));
//			}
//			else {
//				stepl = sqrt(pow((Lpos[s].x() - Lpos[s-1].x()),2) + pow((Lpos[s].y() - Lpos[s-1].y()),2) + pow((Lpos[s].z() - Lpos[s-1].z()),2));
//			}
			
			double Edep_B = BirksAttenuation(Edep[s],Dx[s],charge[s],birks_constant);
			
			// cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
			//	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
			
			if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with
			
			etotL = etotL + Edep_B/2 * exp(-dLeft/att_length) * light_coll;
			etotR = etotR + Edep_B/2 * exp(-dRight/att_length) * light_coll;
			
			etotB = etotB + Edep[s]/2 * exp(-dLeft/att_length) * light_coll;
			etotF = etotF + Edep[s]/2 * exp(-dRight/att_length) * light_coll;
			
			//  cout << "step: " << s << " etotL, etotR, etotB, etotF: " << etotL << ", " << etotR << ", " << etotB << ", " << etotF << endl;
			
			// timeL= timeL + (times[s] + dLeft/veff) / nsteps;
			// timeR= timeR + (times[s] + dRight/veff) / nsteps;
			
			timeL= timeL + (times[s] + dLeft/veff) / nsteps;
			timeR= timeR + (times[s] + dRight/veff) / nsteps;
			
			timeB= timeB + (times[s] + dLeft/veff) / nsteps;
			timeF= timeF + (times[s] + dRight/veff) / nsteps;
			
			if(etotL > 0.) {
				if(s==0 || (time_min[0]>(times[s]+dLeft/veff))) time_min[0]=times[s]+dLeft/veff;
			}
		//.q	      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
			if(etotR > 0.) {
				if(s==0 || (time_min[1]>(times[s]+dRight/veff))) time_min[1]=times[s]+ dRight/veff;
			}
		
			if(etotB > 0.) {
				if(s==0 || (time_min[2]>(times[s]+dLeft/veff))) time_min[2]=times[s]+ dLeft/veff;
			}
			if(etotF > 0.) {
				if(s==0 || (time_min[3]>(times[s]+dRight/veff))) time_min[3]=times[s]+ dRight/veff;
			}
			
		}
	 
		double peL=G4Poisson(etotL*light_yield*sensor_qe);
		double peR=G4Poisson(etotR*light_yield*sensor_qe);
		double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
		double sigmaTR=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peR+1.));
	//	  double sigmaTL=0;
	//	  double sigmaTR=0;
		TDCL=(int) ((time_min[0]+G4RandGauss::shoot(0.,sigmaTL)) * tdc_conv);
		TDCR=(int) ((time_min[1]+G4RandGauss::shoot(0.,sigmaTR)) * tdc_conv);
		if(TDCL<0) TDCL=0;
		if(TDCR<0) TDCR=0;
		ADCL=(int) (peL*sensor_gain*adc_conv + adc_ped);
		ADCR=(int) (peR*sensor_gain*adc_conv + adc_ped);
	//	  cout << "ADCL: " << ADCL << " " << peL << " " << sensor_gain << " " << adc_conv << endl;
		
		
		double peB=G4Poisson(etotB*light_yield*sensor_qe);
		double peF=G4Poisson(etotF*light_yield*sensor_qe);
		double sigmaTB=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peB+1.));
		double sigmaTF=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peF+1.));
	//	  double sigmaTB=0;
	//	  double sigmaTF=0;
		TDCB=(int) ((time_min[2] + G4RandGauss::shoot(0.,sigmaTB)) * tdc_conv);
		TDCF=(int) ((time_min[3] + G4RandGauss::shoot(0.,sigmaTF)) * tdc_conv);
		if(TDCB<0) TDCB=0;
		if(TDCF<0) TDCF=0;
		ADCB=(int) (G4Poisson(etotB*light_yield*sensor_qe)*sensor_gain*adc_conv + adc_ped);
		ADCF=(int) (G4Poisson(etotF*light_yield*sensor_qe)*sensor_gain*adc_conv + adc_ped);

	  //cout << "energy right: " << ADCR / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADCL / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
	  //cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
	  
	  //cout << " Light collection: " << light_coll << endl;
	  
	}
	// closes (Etot > 0) loop
	
	
	if(verbosity>4)
	{
	  cout <<  log_msg << " layer: " << layer    << ", paddle: " << paddle ;
	  cout <<  log_msg << " Etot=" << Etot/MeV << endl;
	  cout <<  log_msg << " TDCL=" << TDCL     << " TDCR=" << TDCR    << " ADCL=" << ADCL << " ADCR=" << ADCR << endl;
	  cout <<  log_msg << " TDCB=" << TDCB     << " TDCF=" << TDCF    << " ADCB=" << ADCB << " ADCF=" << ADCF << endl;
	}
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["layer"]  = layer;
	dgtz["paddle"] = paddle;
	dgtz["adcl"]   = ADCL;
	dgtz["adcr"]   = ADCR;
	dgtz["tdcl"]   = TDCL;
	dgtz["tdcr"]   = TDCR;
	dgtz["adcb"]   = ADCB;
	dgtz["adcf"]   = ADCF;
	dgtz["tdcb"]   = TDCB;
	dgtz["tdcf"]   = TDCF;
	
	return dgtz;
}


vector<identifier>  cormo_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> cormo_HitProcess :: electronicNoise()
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


double cormo_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks)
{
	//Example of Birk attenuation law in organic scintillators.
	//adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
	//
	// Taken from GEANT4 examples advanced/amsEcal and extended/electromagnetic/TestEm3
	//
	double response = destep;
	if (birks*destep*stepl*charge != 0.)
	{
		response = destep/(1. + birks*destep/stepl);
	}
	return response;
}


double cormo_HitProcess::BirksAttenuation2(double destep,double stepl,int charge,double birks)
{
	//Extension of Birk attenuation law proposed by Chou
	// see G.V. O'Rielly et al. Nucl. Instr and Meth A368(1996)745
	// 
	//
	double C=9.59*1E-4*mm*mm/MeV/MeV;
	double response = destep;
	if (birks*destep*stepl*charge != 0.)
	{
		response = destep/(1. + birks*destep/stepl + C*pow(destep/stepl,2.));
	}
	return response;
}


map< string, vector <int> >  cormo_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > cormo_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double cormo_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}



