// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "cnd_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> cnd_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	string hd_msg = " > cnd hit process";
	
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);
	
	int layer  = identity[0].id;
	int paddle = identity[1].id;
	
	// Energy propagation:
	
	double Lg[3];                                 // lengths along light-guides to PMT (from design drawings)
	Lg[0] = 139.56*cm;                            // inner layer
	Lg[1] = 138.89*cm;
	Lg[2] = 136.88*cm;                            // outer layer
	
	double att_length = 1.5*m;                    // light attenuation length in scintillator
	double att_length_lg = 9.5*m;                 // light attenuation length in long light-guides
	
	double sensor_surface = pow(2.5*cm,2)*pi;     // X-sectional area of PMT, assume radius of 2.5 cm.
	double paddle_xsec = 0.;                       // cross-sectional area of the paddle
	if (layer == 1) paddle_xsec = 22.8*cm*cm;     // from the geometry files
	else if (layer == 2) paddle_xsec = 26.4*cm*cm;
	else if (layer == 3) paddle_xsec = 27.7*cm*cm;
	double light_coll;                            // ratio of photo_sensor area over paddle section, times optical coupling ~ light collection efficiency
	if (sensor_surface < paddle_xsec) light_coll = 0.7 * sensor_surface / paddle_xsec;
	else light_coll = 0.7;                        // to make sure sensor_surface / paddle_xsec doesn't go over 1.
	
	double uturn[3];                               // fraction of energy which makes it through the u-turn light-guides (based on cosmic tests)
	uturn[0] = 0.65;                               // inner layer
	uturn[1] = 0.6;
	uturn[2] = 0.5;                                // outer layer
	
	double light_yield = 10000/MeV;               // number of optical photons produced in the scintillator per MeV of deposited energy
	double sensor_qe = 0.2;                       // photo sensor quantum efficiency
	double sensor_gain = 0.24;                    // gain of the photo sensor in pC/(#p.e.); it defines the conversion from photoelectrons to charge:
						      // for a pmt gain of 1.5*10^6, this factor is equal to 1.5*10^6*1.6*10^-19 C = 0.24 pC/(#p.e.)
	double signal_split = 0.5;                    // signal is split into two, going to QDC and discriminators.
	double adc_conv = 10.;                        // conversion factor from pC to ADC (typical sensitivy of CAEN VME QDC is of 0.1 pC/ch)
	double adc_ped = 3.;                          // ADC Pedestal
	
	
	// Time of signal:
	
	double veff = 16*cm/ns;                       // light velocity in scintillator
	
	double t_u[3];                                // time it takes for light to travel round u-turn lightguide for the three layers (based on cosmic tests)
	t_u[0] = 0.6*ns;                              // inner layer
	t_u[1] = 1.2*ns;
	t_u[2] = 1.7*ns;                              // outer layer
	
	double sigmaTD = 0.14*ns/sqrt(MeV);           // time smearing factor (estimated from tests at Orsay), same paddle as hit (in ns/sqrt(MeV)).
	double sigmaTN = 0.14*ns/sqrt(MeV);           // time smearing factor, neighbouring paddle to the hit one, units as above.
	
	double tdc_conv = 40/ns;                      // TDC conversion factor (1/0.025ns), channels per per ns.
	
	
	// Get the paddle length: in CND paddles are along z
	double length = aHit->GetDetector().dimensions[0]; // this is actually the half-length!
	
	// Get info about detector material to eveluate Birks effect
	double birks_constant=aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
	
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<double>      Edep = aHit->GetEdep();
	// Charge for each step
	vector<int> charge = aHit->GetCharges();
	vector<double> times = aHit->GetTime();
	vector<double> dx = aHit->GetDx();
	
	unsigned nsteps = times.size();
	
	
	// Variables for each step:
	
	double dDir = 0.;      // distance travelled by the light along the paddle UPSTREAM (direct light)
	double dNeigh = 0.;    // distance travelled by the light along the paddle DOWNSTREAM (and then through the neighbour paddle)
	double Edep_B = 0.;    // Energy deposit, scaled in accordance with Birk's effect
	
	double e_dir = 0.;     // attenuated energy as it arrives at the two PMTs (one coupled to the hit paddle, one to its neighbour)
	double e_neigh = 0.;
	
	
	// Variables for each hit:
	
	//  int flag_counted = 0;   // Flag to specify that the time of hit has already been determined for this hit (1: yes, 0: no)
	
	double et_D = 0.;      // variables to hold total energy*time values for each step, for the two PMTs
	double et_N = 0.;
	
	double etotD = 0.;  // total energy of hit propagated to the PMT connected to the hit paddle
	double etotN = 0.;  // total energy of hit propagated to downstream end of the hit paddle, round u-turn, along neighbouring paddle and light-guide and into PMT
	double timeD = 0.;  // hit times measured at the upstream edges of the two paddles
	double timeN = 0.;
	
	int ADCD = 0;
	int ADCN = 0;
	int TDCD = 4096;    // max value of the ADC readout
	int TDCN = 4096;
	
	
	if(tInfos.eTot>0)
	{
		for(unsigned int s=0; s<nsteps; s++)
		{
	  
			// Distances travelled through the paddles to the upstream edges (of the hit paddle and of its coupled neighbour):
			dDir = length + Lpos[s].z();
			dNeigh = 3*length - Lpos[s].z();
	  
			//  cout << "\n Distances: " << endl;
			//  cout << "\t dDir, dNeigh: " << dDir << ", " << dNeigh << ", " << endl;
	  
			// apply Birks effect
			Edep_B = BirksAttenuation(Edep[s],dx[s],charge[s],birks_constant);
	  
			//cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " Dx=" << dx[s]
			//     << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
	  
			// Calculate attenuated energy which will reach both PMTs:
			e_dir   = (Edep_B/2) * exp(-dDir/att_length - Lg[layer-1]/att_length_lg) * light_coll;
			e_neigh = (Edep_B/2) * exp(-dNeigh/att_length - Lg[layer-1]/att_length_lg) * uturn[layer-1] * light_coll;
	  
			// Integrate energy over entire hit. These values are the output and will be digitised:
			etotD = etotD + e_dir;
			etotN = etotN + e_neigh;
	  
			//  cout << "step: " << s << " etotD, etotN: " << etotD << ", " << etotN << endl;
	  
	  
			/****** Time of hit calculation *******/
			// In all the methods below, the assumption is that the time taken to travel along the light-guides,
			// through PMTs and the electronics to the ADC units is the same for all paddles, therefore this time offset is ignored.
	  
			// This takes average time of all the steps with energy deposit above 0:
			//if (Edep[s] > 0.){
			//	timeD = timeD + (times[s] + dDir/veff) / nsteps;
			//	timeN = timeN + (times[s] + dNeigh/veff + t_u[layer-1]) / nsteps;
			//}
	  
			// This takes the time of the first step (in order of creation) with energy deposit above 0:
			// if (flag_counted == 0 && Edep[s] > 0.){
			//	timeD = times[s] + dDir/veff;
			//	timeN = times[s] + dNeigh/veff + t_u[layer-1];
			//	flag_counted = 1;                               // so that subsequent steps are not counted in the hit
			// }
	  
			// This calculates the total energy * time value at edges of both paddles,
			// will be used to get the energy-weighted average time of hit (should correspond roughly to the peak energy deposit):
			et_D = et_D + ((times[s] + dDir/veff) * e_dir);
			et_N = et_N + ((times[s] + dNeigh/veff + t_u[layer-1]) * e_neigh);
	  
		}   // close loop over steps s
		
		
		
		/**** The following calculates the time based on energy-weighted average of all step times ****/
		
		timeD = et_D / etotD;      // sum(energy*time) /  sum(energy)
		timeN = et_N / etotN;
		
		
		/******** end timing determination ***********/
		
		// cout << "Reconstructed time (ns): " << (timeD + timeN)/2. - 2.*length/veff - t_u[layer-1]/2. << endl;
		// cout << "Reconstructed z (cm, wrt paddle center): " << (length + t_u[layer-1]*veff/2. - (timeN - timeD)*veff/2.)/10. << endl;
		// cout << "etotD, etotN (in MeV): " << etotD << ", " << etotN << endl;
		// cout << "timeD, timeN (in ns): " << timeD << ", " << timeN << endl;
		
		// Check the full output:
		// cout << "Total steps in this hit: " << nsteps << endl;
		// for (int s=0; s<nsteps; s++)
		// 	{
		// 	  cout << "\n Edep (in MeV): "  << Edep[s] << " time (in ns): " << times[s] << " Lpos-x (in mm): "
		// 	       << Lpos[s].x() << " Lpos-y: " << Lpos[s].y() << " Lpos-z: " << Lpos[s].z() << " pos-x (in mm): "
		// 	       << pos[s].x() << " pos-y: "  << pos[s].y()  << " pos-z: " << pos[s].z() << " etotD (in MeV): "
		// 	       << etotD << " etotN: " << etotN  << " timeD (in ns): " << timeD << " timeN: " << timeN << endl;
		//   }
		
		
		/**** Actual digitisation happens here! *****/
		
		if (etotD > 0.)
		{
			TDCD = (int) ((timeD + G4RandGauss::shoot(0.,sigmaTD/sqrt(etotD))) * tdc_conv);
			ADCD = (int) (G4Poisson(etotD*light_yield*sensor_qe)*signal_split*sensor_gain*adc_conv + adc_ped);
		}
		if (etotN > 0.)
		{
			TDCN = (int) ((timeN + G4RandGauss::shoot(0.,sigmaTN/sqrt(etotN))) * tdc_conv);
			ADCN = (int) (G4Poisson(etotN*light_yield*sensor_qe)*signal_split*sensor_gain*adc_conv + adc_ped);
		}
		
		if(TDCD < 0) TDCD = 0;
		else if (TDCD > 4096) TDCD = 4096;
		if(TDCN < 0) TDCN = 0;
		else if(TDCN > 4096) TDCN = 4096;
		
		if(ADCD < 0) ADCD = 0;
		if(ADCN < 0) ADCN = 0;
		
	}  // closes tInfos.eTot>0
	
	
	if(verbosity>4)
	{
		cout <<  hd_msg << " layer: " << layer    << ", paddle: " << paddle  << " x=" << tInfos.x/cm << "cm, y=" << tInfos.y/cm << "cm, z=" << tInfos.z/cm << "cm" << endl;
		cout <<  hd_msg << " Etot=" << tInfos.eTot/MeV     << "MeV, average time=" << tInfos.time  << "ns"  << endl;
		cout <<  hd_msg << " TDCD= " << TDCD     << ", TDCN= " << TDCN    << ", ADCD= " << ADCD << ", ADCN= " << ADCN << endl;
	}
	
	
	dgtz["hitn"]   = hitn;
	dgtz["layer"]  = layer;
	dgtz["paddle"] = paddle;
	dgtz["ADCD"]   = ADCD;
	dgtz["ADCN"]   = ADCN;
	dgtz["TDCD"]   = TDCD;
	dgtz["TDCN"]   = TDCN;
	
	return dgtz;
}


vector<identifier>  cnd_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


double cnd_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks)
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


map< string, vector <int> >  cnd_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}





