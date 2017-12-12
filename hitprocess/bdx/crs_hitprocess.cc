// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "crs_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

map<string, double> crs_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;
	vector<identifier> identity = aHit->GetId();

	int sector = identity[0].id; 
	int xch  = identity[1].id;
	int ych = identity[2].id; 
		
	// Digitization Parameters
//	double light_yield=1600/MeV;             // number of optical photons pruced in the scintillator per MeV of deposited energy
//	double att_length=108*cm;                 // light at tenuation length
//	double sensor_surface=pow(3.5*cm,2)*pi;   // area of photo sensor
//	double paddle_surface=pow(10*cm,2);       // paddle surface
//	double light_coll=sensor_surface/paddle_surface;  // ratio of photo_sensor area over paddle section ~ light collection efficiency

//	double sensor_qe=0.15;                     // photo sensor quantum efficiency
//	double sensor_gain=0.472*1.6*1.275;         // pmt gain x electron charge in pC (4x10^6)x(1.6x10^-7) ->val * 1.6 * 0.6 = 0,963   modifica 21/01/10
//	double adc_conv=10;                       // conversion factor from pC to ADC (typical sensitivy of CAEN VME QDC is of 0.1 pC/ch)
//	double adc_ped=0;                         // ADC Pedestal
//	double veff=13*cm/ns;                     // light velocity in scintillator
//	double tdc_conv=1/0.001/ns;               // TDC conversion factor
	
	
	// initialize ADC and TDC
//	double etotL = 0;
//	double etotR = 0;
//	double timeL = 0;
//	double timeR = 0;
//	int ADCL = 0;
//	int ADCR = 0;
//	int TDCL = 4096;
//	int TDCR = 4096;
	
//	double etotB = 0; // signal propagating to backward end of paddle hit happened in
//	double etotF = 0; // signal propagating to forward end of the paddle the hit happened in, round u-turn and along neighbouring paddle
//	double timeB = 0;
//	double timeF = 0;
//	int ADCB = 0;
//	int ADCF = 0;
 
 
    // Parameter for BaBar Crystal
    double integration_frac=1.; // integration time SEE DIGI ROUTINE = 1us 67% of the all signal (0.8us->58% 1us->67% 2us->92% 3us->100%)
    double optical_coupling=0.8;
    double light_yield_crs=50000*integration_frac/MeV;
    double att_length_crs=600*cm;
    double sensor_surface_crs=pow(0.3*cm,2);
    // ! requires to be matched with the Babar crystal individual geometry (short side ave = (4.1+4.7)/2=4.4)
    double redout_surface_crs=pow(4.7*cm,2);
    double light_coll_crs=sensor_surface_crs/redout_surface_crs;
    double sensor_qe_crs=0.22; // 24% sipm 100um 22% sipm 25um 35% sipm 50um
//    double sensor_pe_crs=20; //20mV*100ns/50Ohm/2 -> 1 pe = 20 pC
//    double sensor_gain_crs=1;
    // ! requires to be matched with the Babar crystal individual geometry (32,5 cm)
    double length_crs;
    double etotL_crs = 0; //L= Large side redout
    double etotR_crs = 0; //R= short side redout
    double timeL_crs = 0;
    double timeR_crs = 0;
    double veff_crs=30/1.8*cm/ns;                     // light velocity in crystal
//    double adc_conv_crs=1;                       // conversion factor from pC to ADC (typical sensitivy of CAEN VME QDC is of 0.1 pC/ch)
//    double adc_ped_crs=0;                         // ADC Pedestal
    double tdc_conv_crs=1./ns;               // TDC conversion factor
    double  T_offset_crs=0*ns;
    double ADCL_crs = 0;
    double ADCR_crs = 0;
    double TDCL_crs = 4096;
    double TDCR_crs = 4096;
    double TDCB = 4096;
    
    
    
	// Get the paddle length: in crs paddles are along y
//	double length = aHit->GetDetector().dimensions[2];
	//double length = 20*cm;
	
	// Get info about detector material to eveluate Birks effect
	double birks_constant=aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
	// cout << " Birks constant is: " << birks_constant << endl;
	// cout << aHit->GetDetector().GetLogical()->GetMaterial()->GetName() << endl;
    // forcing Birks for CsI(Tl) from arxiv.org/pdf/0911.3041 and checked with alpha
    birks_constant=3.2e-3;
	
//	double time_min[4] = {0,0,0,0};
    double time_min_crs[4] = {0,0,0,0};
	
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double>      Edep = aHit->GetEdep();
	vector<G4double>      Dx   = aHit->GetDx();
    length_crs = aHit->GetDetector().dimensions[4];
    //cout<<length_crs<< endl;
    
	// Charge for each step
	vector<int> charge = aHit->GetCharges();
	vector<G4double> times = aHit->GetTime();
    //vector<string> theseMats = aHit->GetMaterials();
	
	unsigned int nsteps = Edep.size();
//	double       Etot   = 0;
    double       Etot_crs   = 0;
    {
        for(unsigned int s=0; s<nsteps; s++)
        {
//        Etot = Etot + Edep[s];
        Etot_crs = Etot_crs + Edep[s];
        }
    }
    
    double Etot_B_crs=0;
	if(Etot_crs>0)
	{
	  for(unsigned int s=0; s<nsteps; s++)
		{
			// Distances from Downstrem (R), Upstream (L)
            // Sipm's on D side (R)
            double dLeft_crs  = length_crs + Lpos[s].z();//Upstream
            double dRight_crs = length_crs - Lpos[s].z();//Downstream
            //cout<<length_crs<<"  "<< Lpos[s].z()<< endl;
            //cout<<dLeft_crs<<"  "<< dRight_crs<< endl;
            
//            cout << " material " << theseMats[s] << endl;
            
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
			
//			double Edep_B = BirksAttenuation(Edep[s],Dx[s],charge[s],birks_constant);
            
            double Edep_B_crs = BirksAttenuation(Edep[s],Dx[s],charge[s],birks_constant);
            Etot_B_crs = Etot_B_crs + Edep_B_crs;

			// cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
			//	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
			
//			if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with
			if (light_coll_crs > 1) light_coll_crs = 1.;
            
//			etotL = etotL + Edep_B/2 * exp(-dLeft/att_length) * light_coll;
//			etotR = etotR + Edep_B/2 * exp(-dRight/att_length) * light_coll;

            etotL_crs = etotL_crs + Edep_B_crs/2 * exp(-dLeft_crs/att_length_crs) * light_coll_crs;
            etotR_crs = etotR_crs + Edep_B_crs/2 * exp(-dRight_crs/att_length_crs) * light_coll_crs;

            
//			etotB = etotB + Edep[s]/2 * exp(-dLeft/att_length) * light_coll;
//			etotF = etotF + Edep[s]/2 * exp(-dRight/att_length) * light_coll;
			
			//  cout << "step: " << s << " etotL, etotR, etotB, etotF: " << etotL << ", " << etotR << ", " << etotB << ", " << etotF << endl;
			
			// timeL= timeL + (times[s] + dLeft/veff) / nsteps;
			// timeR= timeR + (times[s] + dRight/veff) / nsteps;
			
//			timeL= timeL + (times[s] + dLeft/veff) / nsteps;
//			timeR= timeR + (times[s] + dRight/veff) / nsteps;

            
            timeL_crs= timeL_crs + (times[s] + dLeft_crs/veff_crs) / nsteps;
            timeR_crs= timeR_crs + (times[s] + dRight_crs/veff_crs) / nsteps;

//			timeB= timeB + (times[s] + dLeft/veff) / nsteps;
//			timeF= timeF + (times[s] + dRight/veff) / nsteps;
			
//			if(etotL > 0.) {
//				if(s==0 || (time_min[0]>(times[s]+dLeft/veff))) time_min[0]=times[s]+dLeft/veff;
//			}
		//.q	      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
//			if(etotR > 0.) {
//				if(s==0 || (time_min[1]>(times[s]+dRight/veff))) time_min[1]=times[s]+ dRight/veff;
//			}
            
            
            if(etotL_crs > 0.) {
                if(s==0 || (time_min_crs[0]>(times[s]+dLeft_crs/veff_crs))) time_min_crs[0]=times[s]+dLeft_crs/veff_crs;
            }
            //.q	      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
            if(etotR_crs > 0.) {
                if(s==0 || (time_min_crs[1]>(times[s]+dRight_crs/veff_crs))) time_min_crs[1]=times[s]+ dRight_crs/veff_crs;
            }
          
            //cout<<time_min_crs[1]<<"  "<< dRight_crs<< endl;
		
//			if(etotB > 0.) {
//				if(s==0 || (time_min[2]>(times[s]+dLeft/veff))) time_min[2]=times[s]+ dLeft/veff;
//			}
//			if(etotF > 0.) {
//				if(s==0 || (time_min[3]>(times[s]+dRight/veff))) time_min[3]=times[s]+ dRight/veff;
//			}
//
		}
	 
        // cout << "Etot " << Etot_crs  << " EtotB " << Etot_B_crs  << " steps " << nsteps<< endl;
//		double peL=G4Poisson(etotL*light_yield*sensor_qe);
//		double peR=G4Poisson(etotR*light_yield*sensor_qe);
//        double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
//        double sigmaTR=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peR+1.));

        // Keeping the L side for future pindiode parametrization
        //double peL_crs=G4Poisson(etotL_crs*light_yield_crs*sensor_qe_crs*optical_coupling);
        //double sigmaTL_crs=sqrt(pow(5.*nanosecond,2.)+pow(10.*nanosecond,2.)/(peL_crs/10.+1.));
        // All time spread added in the function WaveForm
        
        // sampling the short side readout and getting timing out of the digitized WF


        // cout << "Etot " << peL_crs   << " EtotB " << rr  << endl;
        
	//	  double sigmaTL=0;
	//	  double sigmaTR=0;
//		TDCL=(int) ((time_min[0]+G4RandGauss::shoot(0.,sigmaTL)) * tdc_conv);
//		TDCR=(int) ((time_min[1]+G4RandGauss::shoot(0.,sigmaTR)) * tdc_conv);
//		if(TDCL<0) TDCL=0;
//		if(TDCR<0) TDCR=0;
//		ADCL=(int) (peL*sensor_gain*adc_conv + adc_ped);
//		ADCR=(int) (peR*sensor_gain*adc_conv + adc_ped);
	//	  cout << "ADCL: " << ADCL << " " << peL << " " << sensor_gain << " " << adc_conv << endl;

        
        
//        TDCL_crs=(int) ((time_min_crs[0]+T_offset_crs+G4RandGauss::shoot(0.,sigmaTL_crs)) * tdc_conv_crs);
//        TDCR_crs=(int) ((time_min_crs[1]+T_offset_crs+G4RandGauss::shoot(0.,sigmaTR_crs)) * tdc_conv_crs);
//        if(TDCL_crs<0) TDCL_crs=0;
//        if(TDCR_crs<0) TDCR_crs=0;
//        ADCL_crs=(int) (peL_crs*sensor_gain_crs*adc_conv_crs + adc_ped_crs);
//        ADCR_crs=(int) (peR_crs*sensor_gain_crs*adc_conv_crs + adc_ped_crs);
        
        // Crystal readout
        //double test=WaveForm(peR_crs,TDCR_crs);
        double* test;
        double tim;
        double peR_int_crs;
        double peR_crs;
        int Nsamp_int=250;  // 1us
        
        //      SIPM 1
        peR_crs=G4Poisson(etotR_crs*light_yield_crs*sensor_qe_crs*optical_coupling);
        test=WaveForm(peR_crs,&tim);
        // integrating over the integration time (each sample 4.ns, see digi routine)
        //for(unsigned int s=0; s<1000; s++)cout << s  << " " <<  test[s] << endl ;
        peR_int_crs=0.;
        for(int s=0; s<Nsamp_int; s++) peR_int_crs=peR_int_crs+test[s];
        //cout << "TDCR: " << tim << " Npe before digi " << peR_crs<< " Npe avter digi " <<peR_int_crs <<endl;
        TDCR_crs=int(tim);
        ADCR_crs=int(peR_int_crs);
        //ADCL_crs=etotR_crs*1000;
        //cout << "TDCR: " << tim << " Npe before digi " << peR_crs<< " Npe avter digi " <<peR_int_crs <<endl;
        
        //Assigning to TDCB the usual timing seen by sipm1 (TDCR)
        double sigmaTR_crs=sqrt(pow(5.*nanosecond,2.)+pow(10.*nanosecond,2.)/(peR_crs/10.+1.));
        sigmaTR_crs=0.;
        TDCB=((time_min_crs[1]+T_offset_crs+G4RandGauss::shoot(0.,sigmaTR_crs)) * tdc_conv_crs);
        //cout <<time_min_crs[1]<<"  "<<TDCB<<endl;
        
        
        //      SIPM 2
        double sensor_qe_crs_sipm2=0.35;//PDE sipm2
        double optical_coupling_sipm2=optical_coupling;// assuming the two optical coupling are the same
        peR_crs=G4Poisson(etotR_crs*light_yield_crs*sensor_qe_crs_sipm2*optical_coupling_sipm2);
        test=WaveForm(peR_crs,&tim);
        // integrating over the integration time (each sample 4.ns, see digi routine)
        //for(unsigned int s=0; s<1000; s++)cout << s  << " " <<  test[s] << endl ;
        peR_int_crs=0.;
        for(int s=0; s<Nsamp_int; s++) peR_int_crs=peR_int_crs+test[s];
        //cout << "TDCR: " << tim << " Npe before digi " << peR_crs<< " Npe avter digi " <<peR_int_crs <<endl;
        TDCL_crs=int(tim); // assigning to L the sipm2
        ADCL_crs=int(peR_int_crs);// assigning to L the sipm2
        //ADCL_crs=etotR_crs*1000;
        //cout << "TDCR: " << tim << " Npe before digi " << peR_crs<< " Npe avter digi " <<peR_int_crs <<endl;
        
       
        
        
//		double peB=G4Poisson(etotB*light_yield*sensor_qe);
//		double peF=G4Poisson(etotF*light_yield*sensor_qe);
//		double sigmaTB=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peB+1.));
//		double sigmaTF=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peF+1.));
	//	  double sigmaTB=0;
	//	  double sigmaTF=0;
//		TDCB=(int) ((time_min[2] + G4RandGauss::shoot(0.,sigmaTB)) * tdc_conv);
//		TDCF=(int) ((time_min[3] + G4RandGauss::shoot(0.,sigmaTF)) * tdc_conv);
//		if(TDCB<0) TDCB=0;
//		if(TDCF<0) TDCF=0;
//		ADCB=(int) (G4Poisson(etotB*light_yield*sensor_qe)*sensor_gain*adc_conv + adc_ped);
//		ADCF=(int) (G4Poisson(etotF*light_yield*sensor_qe)*sensor_gain*adc_conv + adc_ped);

	  //cout << "energy right: " << ADCR / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADCL / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
	  //cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
	  
	  //cout << " Light collection: " << light_coll << endl;
	  
	}
	// closes (Etot > 0) loop
	
	if(verbosity>4)
	{
	  cout <<  log_msg << " xch: " << xch    << ", ych: " << ych ;
	  cout <<  log_msg << " Etot=" << Etot_crs/MeV << endl;
	  cout <<  log_msg << " TDCL=" << TDCL_crs     << " TDCR=" << TDCR_crs    << " ADCL=" << ADCL_crs << " ADCR=" << ADCR_crs << endl;
	  //cout <<  log_msg << " TDCB=" << TDCB     << " TDCF=" << TDCF    << " ADCB=" << ADCB << " ADCF=" << ADCF << endl;
	}
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["xch"]  = xch;
	dgtz["ych"] = ych;
	dgtz["adcl"]   = ADCL_crs;//Downstream side large cell sipm
	dgtz["adcr"]   = ADCR_crs;//Downstream side small cell sipm
	dgtz["tdcl"]   = TDCL_crs;//Downstream side large cell sipm
	dgtz["tdcr"]   = TDCR_crs;//Downstream side small cell sipm
	dgtz["adcb"]   = 0;
	dgtz["adcf"]   = 0;
	dgtz["tdcb"]   = TDCB*1000.;//original time in ps
	dgtz["tdcf"]   = 0;
	
	return dgtz;
}


vector<identifier>  crs_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}



double crs_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks)
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


double crs_HitProcess::BirksAttenuation2(double destep,double stepl,int charge,double birks)
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


map< string, vector <int> >  crs_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}




//double crs_HitProcess::WaveForm(double npe, double time)
double* crs_HitProcess::WaveForm(double npe, double* time)
{
    double c=exp(-2.);
//    double Time;
    double t; // time in usec
    double WF;
    double y;
    double rr;
    int it;
    int Nch_digi=800; //Number of cjannel for the digitizer
    static double WFsample[1000]; //Needs to be >  Nch_digi+size of the response to the single pe
    double smp_t=4./1000. ;// Assuming fADC sampling at 250 MHz 1sample every 4ns

    // double p[6] = {0.14,-3.5,2.5,-2.,0.5,-1.2};
    double p[6] = {0.,0.680,0.64,3.34,0.36,0.};// Babar CsI paprameters: fast component(in us), % fast, slow comp(in us), % slow
    // double p1[6] = {0.33,-0.04,3.45,-0.05,2.5,-0.045};

    double tau=15.; // ampli response time constant (in ns)
    double t0=0.01; // t0 starting time (in ns)
    double area=(tau/c/2.);
    double A=1./area; // amplitude at mnax (55.41 to have it normalized to integral=1, otherwise the max is at 1)
//    double threshold=10.*1./area/smp_t/1000.; //time threshold in pe - 1/55.41/smp_t*1000. is the funct max -

    
    double t_spread= 1.*0.000;// pream time spread in us
    double A_spread= 1.*0.05*A;// pream amp spread (in fraction of 1pe amplitude = A)
    double func=0.;
    // Building the waveform
    for(unsigned int s=0; s<1000; s++)
    {
       WFsample[s]=0;
    }
    // Building the response to a single pe (preamps response)
    double AmpWF[80];
    for(unsigned int s=0; s<80; s++)
    {t=1000.*s*smp_t;
        // parametrization of preamp out time is in ns (rise ~10ns decay~80ns) sampled in 160ns or 40 samples
        //func=1./411.5*((1-exp(p1[0]+p1[1]*t))*exp(p1[2]+p1[3]*t)+exp(p1[4]+p1[5]*t)));
        func=(t-t0)*(t-t0)*exp(-(t-t0)/tau)*A/(4*tau*tau*c)*0.5*(abs(t-t0)/(t-t0)+1);
        // spreading amplitude by apli noise
         AmpWF[s]=smp_t*1000.*func;
    }
    
    // fraction of pe in Nch_digi
    double frac=1-((p[2]*exp(-smp_t*Nch_digi/p[1])+p[4]*exp(-smp_t*Nch_digi/p[3])));
    int Npe = frac*npe;
//    cout << "Npe  " << Npe<<"  "<<npe << endl;
    
    //Npe=1.;
    //Npe=Npe/100.;
    for(int s=1; s<=Npe; s++){
     y=1.;
     WF=0.;
     while (y > WF) {
       rr=(rand() % 1000000+1)/1000000.; // rnd number between 0-1
       t= Nch_digi*smp_t*rr; // extracting over 5000 samples range (5000x4ns=20us)
       //WF= 1./5.15*((1-exp(p[0]+p[1]*t))*exp(p[2]+p[3]*t)+exp(p[4]+p[5]*t));
       WF= (p[2]/p[1]*exp(-t/p[1])+p[4]/p[3]*exp(-t/p[3]))/(p[2]/p[1]+p[4]/p[3]);
       rr=(rand() % 10000000+1)/10000000.; // rnd number between 0-1
       y=rr;
       //  cout << "WF " << WF   << " rnd " << y  << endl;
     }
     // spreading time and amplitude of the ampli signal
//        it=t/smp_t;
//        cout << "t " << t << "it " << it        << endl;

     t=G4RandGauss::shoot(t,t_spread);
     if (t<0.) t=0.;
     it=t/smp_t;
//      cout << "t " << t << "it " << it        << endl;
     for(int s=0; s<80; s++)
      { t=1000.*s*smp_t;
          func=AmpWF[s];
          func=G4RandGauss::shoot(func,A_spread);
          if((s+it)<Nch_digi) WFsample[s+it]=WFsample[s+it]+ func;
      }
     }
    
    //for(unsigned int s=0; s<200; s++) cout << s  << " " <<  WFsample[s] << endl ;

        // mimicking a CF discriminatorm at 1/3 of the max signal
    *time=0.;
    double time_max=-100;
    int s=0;
    int s_time_max=0;
    while (time_max<WFsample[s]){
        time_max=1/2.*(WFsample[s+1]+WFsample[s]);
        s_time_max=s;
        *time=1000.*smp_t*s_time_max/3.;
        s++;
    }
   // cout<<s_time_max<<"  "<< time_max<< "  "<<*time <<endl;
    
/* // mimicking a FixedT discriminatorm
     for(unsigned int s=0; s<1000; s++)
     {
      //cout << s  << " " <<  WFsample[s] << endl ;
         //look for the max
         
         if(WFsample[s]>threshold)
         {*time=1000.*s*smp_t; //time in ns
             break;
         }
    }
 */
    
    // get timing from digitized WF
    
    // cout <<  "routine " <<  WFsample[50] << endl ;
    //cout << "Time  " << threshold<< "sample  "<<s<< "time  "<< *time << endl;

    return WFsample;

}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> crs_HitProcess :: electronicNoise()
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







// - charge: returns charge/time digitized information / step
map< int, vector <double> > crs_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double crs_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}






