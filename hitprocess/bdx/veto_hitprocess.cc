// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "veto_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> veto_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector  = identity[0].id;
	int veto_id = identity[1].id;
	int channel = identity[2].id;
		
	// Digitization Parameters
	double light_yield=10000/MeV;             // number of optical photons pruced in the scintillator per MeV of deposited energy
	double att_length=200*cm;                 // light at tenuation length
	double sensor_surface=pow(2.5*cm,2)*pi;   // area of photo sensor
	double paddle_surface=40.*cm*2*cm;        // paddle surface
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
	int ADC1 = 0;
	int ADC2 = 0;
	int TDC1 = 4096;
	int TDC2 = 4096;
	
	
	// Get the paddle length: in veto paddles are along y
	double length = aHit->GetDetector().dimensions[0];
//	double length = 50*cm;
	
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
			double dLeft  = length + Lpos[s].x();
			double dRight = length - Lpos[s].x();
			dLeft  = 0;
			dRight = 0;
			
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
			Edep_B=Edep[s];
			
			// cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
			//	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
			
			if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with
			
			etotL = etotL + Edep_B/2 * exp(-dLeft/att_length) * light_coll;
			etotR = etotR + Edep_B/2 * exp(-dRight/att_length) * light_coll;
			
			//  cout << "step: " << s << " etotL, etotR, etotB, etotF: " << etotL << ", " << etotR << ", " << etotB << ", " << etotF << endl;
			
			timeL= timeL + (times[s] + dLeft/veff) / nsteps;
			timeR= timeR + (times[s] + dRight/veff) / nsteps;
			
			if(etotL > 0.) {
				if(s==0 || (time_min[0]>(times[s]+dLeft/veff))) time_min[0]=times[s]+dLeft/veff;
			}
		//.q	      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
			if(etotR > 0.) {
				if(s==0 || (time_min[1]>(times[s]+dRight/veff))) time_min[1]=times[s]+ dRight/veff;
			}
		}
	 
		double peL=G4Poisson(etotL*light_yield*sensor_qe);
		double peR=G4Poisson(etotR*light_yield*sensor_qe);
		peL=(etotL*light_yield*sensor_qe);
		peR=(etotR*light_yield*sensor_qe);
		double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
		double sigmaTR=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peR+1.));
		sigmaTL=0;
		sigmaTR=0;
		TDC1=(int) ((time_min[0]+G4RandGauss::shoot(0.,sigmaTL)) * tdc_conv);
		TDC2=(int) ((time_min[1]+G4RandGauss::shoot(0.,sigmaTR)) * tdc_conv);
		if(TDC1<0) TDC1=0;
		if(TDC2<0) TDC2=0;
		ADC1=(int) (peL*sensor_gain*adc_conv + adc_ped);
		ADC2=(int) (peR*sensor_gain*adc_conv + adc_ped);
		
		if(!(veto_id==2 && (channel==1&&channel==2&&channel==4&&channel==5))) {
			ADC2=0;
			TDC2=4096;
		}
	//	  cout << "ADC1: " << ADC1 << " " << peL << " " << sensor_gain << " " << adc_conv << endl;


	  //cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
	  //cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
	  
	  //cout << " Light collection: " << light_coll << endl;
	  
	}
	// closes (Etot > 0) loop
	
	
	if(verbosity>4)
	{
	  cout <<  log_msg << " veto: " << veto_id   << ", channel: " << channel ;
	  cout <<  log_msg << " Etot=" << Etot/MeV << endl;
	  cout <<  log_msg << " TDC1=" << TDC1     << " TDC2=" << TDC2    << " ADC1=" << ADC1 << " ADC2=" << ADC2 << endl;
	}
	
	dgtz["hitn"]    = hitn;
	dgtz["sector"]  = sector;
	dgtz["veto"]    = veto_id;
	dgtz["channel"] = channel;
	dgtz["adc1"]    = ADC1;
	dgtz["adc2"]    = ADC2;
	dgtz["tdc1"]    = TDC1;
	dgtz["tdc2"]    = TDC2;
	
	return dgtz;
}


vector<identifier>  veto_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}



double veto_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks)
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


double veto_HitProcess::BirksAttenuation2(double destep,double stepl,int charge,double birks)
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


map< string, vector <int> >  veto_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}





