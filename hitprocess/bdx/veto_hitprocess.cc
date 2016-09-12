
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

    

//	double adc_conv=10;                 // conversion factor from pC to ADC (typical sensitivy of CAEN VME QDC is of 0.1 pC/ch)
//	double adc_ped=0;                   // ADC Pedestal
//	double tdc_conv=1*ns;               // TDC conversion factor
	
	
	// initialize ADC and TDC
	double etotL = 0;
	double etotR = 0;
	double timeL = 0;
	double timeR = 0;
    int ADC1 = 0;
    int ADC2 = 0;
    int ADC3 = 0;
    int ADC4 = 0;
    int TDC1 = 4096;
    int TDC2 = 4096;
    int TDC3 = 4096;
    int TDC4 = 4096;
 
    double length = 0;
    // From measurement
    // spe = 0.36pC
    // Cosmic = 84pC (~235pe)
    // Attenuation at 80cm: ~0.85 -> effective att lenght ~350cm
    //Plastic
    double paddle_surface = 0;        // paddle surface
    double light_yield = 0;  // number of optical photons pruced in the scintillator per MeV of deposited energy
    double att_length = 0;               // light at tenuation length
    double veff = 0;            // light velocity in scintillator
     //PMT
    double sensor_surface = 0;   // area of photo sensor
    double sensor_effective_area = 0; // considering only a fraction of the photocathod
    double sensor_qe = 0;                     // photo sensor quantum efficiency
    double sensor_gain = 0;         // pmt gain x electron charge in pC (2.2x10^6)x(1.6x10^-7) -> ~0.36pC or 1 to have pe
    double light_coll = 0; // ratio of photo_sensor area over paddle section ~ light collection efficiency
    double light_guide_att = 0;
    double tL = 0;
    double tR = 0;
    double peL=0.;
    double peR = 0;
    double etot_g4=0.;

    // Proposal
    // sector: run over channels, from 0 to N
    // Oveto == 5
    // channel: run over the position (1=T 2=B 3=R 4=L 5=D 6=U)
    if(veto_id==5)
    {
        {
            double optical_coupling[13]= {0., 0.94,0.57, 0.35, 0.7, 0.094, 0.177, 0.52, 0.75, 0.52, 0.52, 0.38, 1.0 };
            for (int s=0; s<13; s++) optical_coupling[s] = optical_coupling[s]*0.68;
            light_yield=9200/MeV;
            veff=13*cm/ns    ;
            sensor_effective_area=0.9;
            sensor_qe=0.25;
            sensor_gain=1.;
            
            // Upper/lower
            if(channel==1 || channel==2)
            {
                // Get the paddle length: in veto paddles are along z
                length = aHit->GetDetector().dimensions[2];
                double s1=aHit->GetDetector().dimensions[0];
                double s2=aHit->GetDetector().dimensions[1];
                paddle_surface = 2*s1*2*s2;
                sensor_surface=pow(2.5*cm,2)*pi; // 2" pmt R-> (2*2.5/2 TBC)
                att_length=350*cm;
                light_guide_att=1.0;
                // cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
            }
            if(channel==5 || channel==6)
            {
                double s1=aHit->GetDetector().dimensions[0];
                double s2=aHit->GetDetector().dimensions[3];
                paddle_surface = 2*s1*2*s2;// surface perpendicular to the pmt position Surf=XxZ
                sensor_surface=pow((2.5/2)*cm,2)*pi; //
                att_length=400*cm; // longer att lenght to take into account the perpendicular readout
                light_guide_att=0.19; // no light guide
                // cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
            }
            if(channel==3 || channel==4)
            {
                // Get the paddle length: in veto paddles are along y
                length = aHit->GetDetector().dimensions[1];
                double s1=aHit->GetDetector().dimensions[0];
                double s2=aHit->GetDetector().dimensions[2];
                sensor_surface=pow(2.5*cm,2)*pi; // 2" pmt R-> (2*2.5/2 TBC)
                paddle_surface = 2*s1*2*s2;
                att_length=350*cm;
                light_guide_att=1.0;
                // cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
            }
            
            light_coll=sensor_surface/paddle_surface;
            if (sensor_surface>paddle_surface) light_coll=1.;   // no more than the PMT size
            light_coll=optical_coupling[1]*light_coll*sensor_effective_area*light_guide_att; // coupling [1] identical for all
            //cout << " light collo " << light_coll     <<endl;
            
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
                    double dLeft  =-10000.;
                    double dRight =-10000.;
                    
                    // Distances from left, right for upper/lower (along z)
                    if(channel==1 || channel==2)
                    {
                        dLeft  = length - Lpos[s].z();
                        dRight = length + Lpos[s].z();
                    }
                    // Distances from left, right for other OV (along y)
                    
                    if(channel==3 || channel==4)
                    {
                        dLeft  = length + Lpos[s].y();
                        dRight = length - Lpos[s].y();
                    }
                    if(channel==5 || channel==6)
                    {
                        dLeft  = Lpos[s].x();
                        dRight = Lpos[s].y();
                        double dCent=sqrt(dLeft*dLeft+dRight*dRight);
                        dRight =dCent;
                        
                    }
                    
                    
                    
                    // cout << "\n Distances: " << endl;
                    // cout << "\t dLeft, dRight " << dLeft <<  ", " << dRight << endl;
                    
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
                    etot_g4=etot_g4+Edep_B;
                    // cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
                    //	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
                    
                    //if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with
                    
                    etotL = etotL + Edep_B/2 * exp(-dLeft/att_length) * light_coll;
                    etotR = etotR + Edep_B/2 * exp(-dRight/att_length) * light_coll;
                    
                    //			  cout << "step: " << s << " etotL, etotR " << etotL << ", " << etotR  << endl;
                    
                    timeL= timeL + (times[s] + dLeft/veff) / nsteps;
                    timeR= timeR + (times[s] + dRight/veff) / nsteps;
                    
                    if(etotL > 0.) {
                        if(s==0 || (time_min[0]>(times[s]+dLeft/veff))) time_min[0]=times[s]+dLeft/veff;
                    }
                    //      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
                    if(etotR > 0.) {
                        if(s==0 || (time_min[1]>(times[s]+dRight/veff))) time_min[1]=times[s]+ dRight/veff;
                    }
                }
                //cout << " etotR " << etotR   <<  " ; etotL" <<  etotL <<endl;
                peL=G4Poisson(etotL*light_yield*sensor_qe);
                peR=G4Poisson(etotR*light_yield*sensor_qe);
                //cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
                //peL=(etotL*light_yield*sensor_qe);
                //peR=(etotR*light_yield*sensor_qe);
                //cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
                
                double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
                double sigmaTR=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peR+1.));
                //sigmaTL=0;
                //sigmaTR=0;
                tL=(time_min[0]+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
                tR=(time_min[1]+G4RandGauss::shoot(0.,sigmaTR))*1000.;// time in ps
                // Digitization for ADC and QDC not used
                //TDC1=(int) (tL * tdc_conv);
                //TDC2=(int) (tR * tdc_conv);
                //if(TDC1<0) TDC1=0;
                //if(TDC2<0) TDC2=0;
                //ADC1=(int) (peL*sensor_gain*adc_conv + adc_ped);
                //ADC2=(int) (peR*sensor_gain*adc_conv + adc_ped);
                
                //	  cout << "ADC1: " << ADC1 << " " << peL << " " << sensor_gain << " " << adc_conv << endl;
                
                
                //cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
                //cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
                
                //cout << " Light collection: " << light_coll << endl;
                
            }
            // closes (Etot > 0) loop
            
            
            
            // ch 1,3               dRight
            // ch 2,4               dLeft
            // ch 7,8,9,10,11,12    dRight
            // ch 5,6               dRight
            // ignore the other side
            
            if(channel==1 ||  channel==2 )
            { if(sector==0)
                {ADC1=peR;
                 TDC1=tR;}
            else if(sector==1)
                {ADC1=peL;
                 TDC1=tL;}
            }
            if(channel==3 || channel ==4)
            {
                ADC1=peR;
                TDC1=tR;
            }
            if(channel==5 || channel==6)
            {
                ADC1=peR;
                TDC1=tR;
            }

            // cout << " ADC1: " << ADC1    <<  " ; TDC1: " <<  TDC1  << " ;  ADC2: "<< etot_g4*1000. << endl;
            if(verbosity>4)
            {
                cout <<  log_msg << " veto: " << veto_id   << ", channel: " << channel << ", sector: " << sector ;
                cout <<  log_msg << " Etot=" << Etot/MeV << endl;
                cout <<  log_msg << " TDC1=" << TDC1     << " TDC2=" << TDC2    << " ADC1=" << ADC1 << " ADC2=" << ADC2 << endl;
            }

            
        }
    }
    else if(veto_id==4)
    // proposal IV
    {
        double veff=13*cm/ns    ;// TO BE CHECKED
        // scintillator sizes
        double sx=aHit->GetDetector().dimensions[0];
        double sy=aHit->GetDetector().dimensions[1];
        double sz=aHit->GetDetector().dimensions[2];
        
//        double time_min[4] = {0,0,0,0};
        
        vector<G4ThreeVector> Lpos = aHit->GetLPos();
        vector<G4double>      Edep = aHit->GetEdep();
        vector<G4double>      Dx   = aHit->GetDx();
        // Charge for each step
        vector<int> charge = aHit->GetCharges();
        vector<G4double> times = aHit->GetTime();
        unsigned int nsteps = Edep.size();
        double       Etot   = 0;
        double  X_hit_ave=0.;
        double  Y_hit_ave=0.;
        double  Z_hit_ave=0.;
        double  T_hit_ave=0.;
        double dLeft  =-10000.;
        
        for(unsigned int s=0; s<nsteps; s++) Etot = Etot + Edep[s];
        if(Etot>0)
        {
            for(unsigned int s=0; s<nsteps; s++)
            {
                double Edep_B=Edep[s];
                etot_g4=etot_g4+Edep_B;
                // average hit position XYZ
                X_hit_ave=X_hit_ave+Lpos[s].x();
                Y_hit_ave=Y_hit_ave+Lpos[s].y();
                Z_hit_ave=Z_hit_ave+Lpos[s].z();
                // average hit time
                T_hit_ave=T_hit_ave+times[s];
                
                //cout << "X " << Lpos[s].x() << " " << "Y " << Lpos[s].y() << " " << "Z " << Lpos[s].z() << " "<< "T " <<times[s] << " " << endl;
                
            }
            X_hit_ave=X_hit_ave/nsteps;
            Y_hit_ave=Y_hit_ave/nsteps;
            Z_hit_ave=Z_hit_ave/nsteps;
            T_hit_ave=T_hit_ave/nsteps;
            dLeft  =sz-Z_hit_ave;
            timeL= dLeft/veff+T_hit_ave;
            
            
            double *pe_sipm;// response for a mip (2..05 MeV energy released in 1cm thick)

            pe_sipm=IVresponseProposal(channel, X_hit_ave, Y_hit_ave, Z_hit_ave,sx,sy,sz);
            
            ADC1=G4Poisson(pe_sipm[0]*etot_g4/2.05) ; // Scaling for more/less energy release)
            ADC2=G4Poisson(pe_sipm[1]*etot_g4/2.05) ; // Scaling for more/less energy release)
            ADC3=G4Poisson(pe_sipm[2]*etot_g4/2.05) ; // Scaling for more/less energy release)
            ADC4=G4Poisson(pe_sipm[3]*etot_g4/2.05) ; // Scaling for more/less energy release)
            double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
            sigmaTL=0.;
            TDC1=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            TDC2=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            TDC3=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            TDC4=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            
             //cout <<  log_msg << " veto: " << veto_id   << ", channel: " << channel << ", sector: " << sector << endl;
              //cout << "X " << X_hit_ave << " " << "Y " << Y_hit_ave << " " << "Z " << Z_hit_ave << " "<< "T " <<T_hit_ave << " " << endl;
              //cout << "sipm1 " << pe_sipm[0] << " " << "sipm2 " << pe_sipm[1] << " " << "sipm3 " << pe_sipm[2] << " " << "sipm4 " << pe_sipm[3] << " " << endl;
            // cout << "dLeft " << dLeft << " " << "timeL " << timeL << " " << endl;
            
            //cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
            //cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
            
            //cout << " Light collection: " << light_coll << endl;
            
        }
        // closes (Etot > 0) loop
        
    }
    
    // End Proposal
    
    
    
    
    // Outer veto
    if(veto_id==2)
    {
        double optical_coupling[13]= {0., 0.94,0.57, 0.35, 0.7, 0.094, 0.177, 0.52, 0.75, 0.52, 0.52, 0.38, 1.0 };
        for (int s=0; s<13; s++) optical_coupling[s] = optical_coupling[s]*0.68;
        light_yield=9200/MeV;
        veff=13*cm/ns    ;
        sensor_effective_area=0.9;
        sensor_qe=0.25;
        sensor_gain=1.;

    // Upper/lower
        if(channel==1 || channel==2 || channel==3 || channel ==4)
        {
            // Get the paddle length: in veto paddles are along z
            length = aHit->GetDetector().dimensions[2];
             double s1=aHit->GetDetector().dimensions[0];
            double s2=aHit->GetDetector().dimensions[1];
            paddle_surface = 2*s1*2*s2;
            sensor_surface=pow(2.5*cm,2)*pi; // 2" pmt R-> (2*2.5/2 TBC)
            att_length=350*cm;
            light_guide_att=1.0;
           // cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
        }
        if(channel==5 || channel==6)
        {
            double s1=aHit->GetDetector().dimensions[0];
            double s2=aHit->GetDetector().dimensions[3];
            paddle_surface = 2*s1*2*s2;// surface perpendicular to the pmt position Surf=XxZ
            sensor_surface=pow((2.5/2)*cm,2)*pi; //
            att_length=400*cm; // longer att lenght to take into account the perpendicular readout
            light_guide_att=0.19; // no light guide
           // cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
        }
        if(channel==7 || channel==8 || channel==9 || channel ==10 || channel ==11 || channel ==12 )
        {
            // Get the paddle length: in veto paddles are along y
            length = aHit->GetDetector().dimensions[1];
            double s1=aHit->GetDetector().dimensions[0];
            double s2=aHit->GetDetector().dimensions[2];
            sensor_surface=pow(2.5*cm,2)*pi; // 2" pmt R-> (2*2.5/2 TBC)
            paddle_surface = 2*s1*2*s2;
            att_length=350*cm;
            light_guide_att=1.0;
           // cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
        }
        
        light_coll=sensor_surface/paddle_surface;
        if (sensor_surface>paddle_surface) light_coll=1.;   // no more than the PMT size
        light_coll=optical_coupling[channel]*light_coll*sensor_effective_area*light_guide_att;             // Including the coupling efficiency and the pc effective area
        //cout << " light collo " << light_coll     <<endl;
	
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
            double dLeft  =-10000.;
            double dRight =-10000.;
            
			// Distances from left, right for upper/lower (along z)
            if(channel==1 || channel==2 || channel==3 || channel ==4)
            {
			dLeft  = length - Lpos[s].z();
			dRight = length + Lpos[s].z();
            }
            // Distances from left, right for other OV (along y)
           
            if(channel==7 || channel==8 || channel==9 || channel ==10 || channel ==11 || channel ==12 )
            {
                dLeft  = length + Lpos[s].y();
                dRight = length - Lpos[s].y();
            }
            if(channel==5 || channel==6)
            {
                dLeft  = Lpos[s].x();
                dRight = Lpos[s].y();
                double dCent=sqrt(dLeft*dLeft+dRight*dRight);
                dRight =dCent;
                
            }
           
            
            
			// cout << "\n Distances: " << endl;
            // cout << "\t dLeft, dRight " << dLeft <<  ", " << dRight << endl;
            
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
            etot_g4=etot_g4+Edep_B;
			// cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
			//	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
			
			//if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with
			
			etotL = etotL + Edep_B/2 * exp(-dLeft/att_length) * light_coll;
			etotR = etotR + Edep_B/2 * exp(-dRight/att_length) * light_coll;
			
//			  cout << "step: " << s << " etotL, etotR " << etotL << ", " << etotR  << endl;
			
			timeL= timeL + (times[s] + dLeft/veff) / nsteps;
			timeR= timeR + (times[s] + dRight/veff) / nsteps;
			
			if(etotL > 0.) {
				if(s==0 || (time_min[0]>(times[s]+dLeft/veff))) time_min[0]=times[s]+dLeft/veff;
			}
			//      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
			if(etotR > 0.) {
				if(s==0 || (time_min[1]>(times[s]+dRight/veff))) time_min[1]=times[s]+ dRight/veff;
			}
		}
        //cout << " etotR " << etotR   <<  " ; etotL" <<  etotL <<endl;
		peL=G4Poisson(etotL*light_yield*sensor_qe);
        peR=G4Poisson(etotR*light_yield*sensor_qe);
        //cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
        //peL=(etotL*light_yield*sensor_qe);
		//peR=(etotR*light_yield*sensor_qe);
        //cout << " per " << peR   <<  " ; pel" <<  peL <<endl;

		double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
		double sigmaTR=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peR+1.));
		//sigmaTL=0;
		//sigmaTR=0;
		tL=(time_min[0]+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
		tR=(time_min[1]+G4RandGauss::shoot(0.,sigmaTR))*1000.;// time in ps
        // Digitization for ADC and QDC not used
        //TDC1=(int) (tL * tdc_conv);
        //TDC2=(int) (tR * tdc_conv);
		//if(TDC1<0) TDC1=0;
		//if(TDC2<0) TDC2=0;
		//ADC1=(int) (peL*sensor_gain*adc_conv + adc_ped);
		//ADC2=(int) (peR*sensor_gain*adc_conv + adc_ped);
		
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
        
        
        // ch 1,3               dRight
        // ch 2,4               dLeft
        // ch 7,8,9,10,11,12    dRight
        // ch 5,6               dRight
        // ignore the other side
        
        if(channel==1 ||  channel==3 )
        {
            ADC1=peR;
            TDC1=tR;
        }
        if(channel==2 || channel ==4)
        {
            ADC1=peL;
            TDC1=tL;
        }
        if(channel==5 || channel==6)
        {
            ADC1=peR;
            TDC1=tR;
        }
        if(channel==7 || channel==8 || channel==9 || channel ==10 || channel ==11 || channel ==12 )
        {
            ADC1=peR;
            TDC1=tR;
            
         }
       // cout << " ADC1: " << ADC1    <<  " ; TDC1: " <<  TDC1  << " ;  ADC2: "<< etot_g4*1000. << endl;
        
    } // end of OV
    
    
    
    
    
    
    
// INNER VETO
   if(veto_id==1)
    {
        double veff=13*cm/ns    ;// TO BE CHECKED
        // scintillator sizes
//        double sx=aHit->GetDetector().dimensions[0];
//        double sy=aHit->GetDetector().dimensions[1];
        double sz=aHit->GetDetector().dimensions[2];
        
//        double time_min[4] = {0,0,0,0};
        
        vector<G4ThreeVector> Lpos = aHit->GetLPos();
        vector<G4double>      Edep = aHit->GetEdep();
        vector<G4double>      Dx   = aHit->GetDx();
        // Charge for each step
        vector<int> charge = aHit->GetCharges();
        vector<G4double> times = aHit->GetTime();
        unsigned int nsteps = Edep.size();
        double       Etot   = 0;
        double  X_hit_ave=0.;
        double  Y_hit_ave=0.;
        double  Z_hit_ave=0.;
        double  T_hit_ave=0.;
        double dLeft  =-10000.;
        
        for(unsigned int s=0; s<nsteps; s++) Etot = Etot + Edep[s];
        if(Etot>0)
        {
            for(unsigned int s=0; s<nsteps; s++)
            {
                double Edep_B=Edep[s];
                etot_g4=etot_g4+Edep_B;
                // average hit position XYZ
                X_hit_ave=X_hit_ave+Lpos[s].x();
                Y_hit_ave=Y_hit_ave+Lpos[s].y();
                Z_hit_ave=Z_hit_ave+Lpos[s].z();
                // average hit time
                T_hit_ave=T_hit_ave+times[s];
                
        //cout << "X " << Lpos[s].x() << " " << "Y " << Lpos[s].y() << " " << "Z " << Lpos[s].z() << " "<< "T " <<times[s] << " " << endl;

            }
            X_hit_ave=X_hit_ave/nsteps;
            Y_hit_ave=Y_hit_ave/nsteps;
            Z_hit_ave=Z_hit_ave/nsteps;
            T_hit_ave=T_hit_ave/nsteps;
            dLeft  =sz-Z_hit_ave;
            timeL= dLeft/veff+T_hit_ave;
            
            
            double *pe_sipm;// response for a mip (2..05 MeV energy released in 1cm thick)
            pe_sipm=IVresponse(channel, X_hit_ave, Y_hit_ave, Z_hit_ave);
            
            ADC1=G4Poisson(pe_sipm[0]*etot_g4/2.05) ; // Scaling for more/less energy release)
            ADC2=G4Poisson(pe_sipm[1]*etot_g4/2.05) ; // Scaling for more/less energy release)
            ADC3=G4Poisson(pe_sipm[2]*etot_g4/2.05) ; // Scaling for more/less energy release)
            ADC4=G4Poisson(pe_sipm[3]*etot_g4/2.05) ; // Scaling for more/less energy release)
            double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
            sigmaTL=0.;
            TDC1=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            TDC2=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            TDC3=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            TDC4=(timeL+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps


         //   cout << "channel " << channel << endl;
         // cout << "X " << X_hit_ave << " " << "Y " << Y_hit_ave << " " << "Z " << Z_hit_ave << " "<< "T " <<T_hit_ave << " " << endl;
          //  cout << "sipm1 " << pe_sipm[0] << " " << "sipm2 " << pe_sipm[1] << " " << "sipm3 " << pe_sipm[2] << " " << "sipm4 " << pe_sipm[3] << " " << endl;
           // cout << "dLeft " << dLeft << " " << "timeL" << timeL << " " << endl;
            
            //cout << "energy right: " << ADC2 / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E left: " << ADC1 / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
            //cout << "energy forw: " << ADCF / (adc_conv*sensor_gain*sensor_qe*light_yield) << " E back: " << ADCB / (adc_conv*sensor_gain*sensor_qe*light_yield) << endl;
            
            //cout << " Light collection: " << light_coll << endl;
            
        }
        // closes (Etot > 0) loop
        
    }// end of IV
    
    
    //starting paddles
    if(veto_id==3)
    {
            double optical_coupling[3]= {0., 1.,0.37 };
        for (int s=0; s<3; s++) optical_coupling[s] = optical_coupling[s]*0.34;

        light_yield=9200/MeV;
        veff=13*cm/ns    ;
        sensor_surface=pow(1.27*cm,2)*pi; // 1" pmt R-> (2.5/2 TBC)
        sensor_effective_area=0.9;
        
        sensor_qe=0.25;
        sensor_gain=1.;


         // Get the paddle length: in veto paddles are along z
         length = aHit->GetDetector().dimensions[2];
          double s1=aHit->GetDetector().dimensions[0];
          double s2=aHit->GetDetector().dimensions[1];
          paddle_surface = 2*s1*2*s2;
          att_length=350*cm;
          light_guide_att=1.;
        //cout << " lenght: " << length    <<  " optical-coupled surface: " <<  paddle_surface   <<endl;
        
        light_coll=sensor_surface/paddle_surface;
        if (sensor_surface>paddle_surface) light_coll=1.;   // no more than the PMT size
        light_coll=optical_coupling[channel]*light_coll*sensor_effective_area*light_guide_att;             // Including the coupling efficiency and the pc effective area
       // cout << " channel,veto: " << channel << " " << veto_id    <<  " optical-couping: " <<  optical_coupling[channel]  <<  " light coll: " <<  light_coll  <<endl;
        //cout << " light collo " << light_coll     <<endl;
        
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
        //for(unsigned int s=0; s<nsteps; s++) cout << "Energy = " << Edep[s]*1000.*1000. << endl;
        
        if(Etot>0)
        {
            for(unsigned int s=0; s<nsteps; s++)
            {
                double dLeft  =-10000.;
                double dRight =-10000.;
                
                // Distances from left, right for upper/lower (along z)
                    dLeft  = length - Lpos[s].z();
                    dRight = length + Lpos[s].z();
                
                
                //cout << "\n Distances: " << endl;
                // cout << "\t dLeft, dRight " << dLeft <<  ", " << dRight << endl;
                
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
                etot_g4=etot_g4+Edep_B;
                // cout << "\t Birks Effect: " << " Edep=" << Edep[s] << " StepL=" << stepl
                //	  << " PID =" << pids[s] << " charge =" << charge[s] << " Edep_B=" << Edep_B << endl;
                
                //if (light_coll > 1) light_coll = 1.;     // To make sure you don't miraculously get more energy than you started with
                
                etotL = etotL + Edep_B/2 * exp(-dLeft/att_length) * light_coll;
                etotR = etotR + Edep_B/2 * exp(-dRight/att_length) * light_coll;
                
                //			  cout << "step: " << s << " etotL, etotR " << etotL << ", " << etotR  << endl;
                
                timeL= timeL + (times[s] + dLeft/veff) / nsteps;
                timeR= timeR + (times[s] + dRight/veff) / nsteps;
                
                if(etotL > 0.) {
                    if(s==0 || (time_min[0]>(times[s]+dLeft/veff))) time_min[0]=times[s]+dLeft/veff;
                }
                //      cout << "min " << time_min[0] << "min " << times[s]+dLeft/veff << endl;
                if(etotR > 0.) {
                    if(s==0 || (time_min[1]>(times[s]+dRight/veff))) time_min[1]=times[s]+ dRight/veff;
                }
            }
            //cout << " etotR " << etotR   <<  " ; etotL" <<  etotL <<endl;
            peL=G4Poisson(etotL*light_yield*sensor_qe);
            peR=G4Poisson(etotR*light_yield*sensor_qe);
            //cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
            //peL=(etotL*light_yield*sensor_qe);
            //peR=(etotR*light_yield*sensor_qe);
            //cout << " per " << peR   <<  " ; pel" <<  peL <<endl;
            
            double sigmaTL=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peL+1.));
            double sigmaTR=sqrt(pow(0.2*nanosecond,2.)+pow(1.*nanosecond,2.)/(peR+1.));
//            sigmaTL=0;
//            sigmaTR=0;
            tL=(time_min[0]+G4RandGauss::shoot(0.,sigmaTL))*1000.;//time in ps
            tR=(time_min[1]+G4RandGauss::shoot(0.,sigmaTR))*1000.;// time in ps
          //  cout << " tL " << tL   <<  " ; timeL" <<  timeL <<endl;
            // Digitization for ADC and QDC not used
            //TDC1=(int) (tL * tdc_conv);
            //TDC2=(int) (tR * tdc_conv);
            //if(TDC1<0) TDC1=0;
            //if(TDC2<0) TDC2=0;
            //ADC1=(int) (peL*sensor_gain*adc_conv + adc_ped);
            //ADC2=(int) (peR*sensor_gain*adc_conv + adc_ped);
            
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
        
        
        // ch 1,3               dRight
        // ch 2,4               dLeft
        // ch 7,8,9,10,11,12    dRight
        // ch 5,6               dRight
        // ignore the other side
        

            ADC1=peR;
            TDC1=tR;
        // cout << " ADC1: " << ADC1    <<  " ; TDC1: " <<  TDC1  << " ;  ADC2: "<< etot_g4*1000. << endl;

     }// end of paddles
    
    
	dgtz["hitn"]    = hitn;
	dgtz["sector"]  = sector;
	dgtz["veto"]    = veto_id;
	dgtz["channel"] = channel;
	dgtz["adc1"]    = ADC1;// output in pe
	dgtz["adc2"]    = ADC2;//deposited energy in keV
    dgtz["adc3"]    = ADC3;// ignore
    dgtz["adc4"]    = ADC4;// ignore
	dgtz["tdc1"]    = TDC1;// output in ps
	dgtz["tdc2"]    = TDC2;// ignore
    dgtz["tdc3"]    = TDC3;// ignore
    dgtz["tdc4"]    = TDC4;// ignore
	
	return dgtz;
}


vector<identifier>  veto_HitProcess :: processID(vector<identifier> id, G4Step *step, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


double* veto_HitProcess::IVresponse(int channel, double xx, double yy,double zz)
{
    
    // Response of the different IV plastic paddles
    // ch
    //
    //
    static double response[4];

    for(unsigned int s=0; s<4; s++)response[s] = 0.;
    
    if (channel==1)//top
    {
        double x=xx/10.;
        double y=(1058/2.-zz)/10.;
        
        double parm[4][8]={
        {1.99627e+01,	1.64910e-01,	-5.83528e-01,	-7.34483e-03,	-1.25062e-03, 	4.43805e-03,	5.63766e-05, 	1.40682e-05},
        {1.86162e+01,	4.36475e-02,	-6.78752e-02,	-5.47887e-03,	-1.60512e-04,	-2.33958e-02,	5.55285e-05,	-5.94424e-05},
        {1.85966e+01,	1.96301e-01,	1.34868e-01,	-7.66131e-04,	-1.61720e-03,	-1.91598e-02,	-1.76198e-06,	-4.72970e-05},
        {9.73394e+00,	1.56111e-01,	3.27558e-01,	2.45041e-03,	-1.31615e-03,	5.82688e-03,	-1.48528e-05,	2.35177e-05}
        
    };

       for(unsigned int s=0; s<4; s++)
       response[s] = parm[s][7]*x*x*y + parm[s][6]*x*y*y + parm[s][5]*x*x + parm[s][4]*y*y + parm[s][3]*x*y + parm[s][2]*x + parm[s][1]*y + parm[s][0];

   //     cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
        //cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;

    }
    else if (channel==2)//bottom
    {// Assuming an overall size of 42.8 cm with 4 bars of
        double x=-(xx-428/2)/10;
        double y=(1058/2.-zz)/10.;
        
        
        for(unsigned int s=0; s<4; s++) response[s] =0.;
        if (x<10)           response[0]= (-0.000303034)*y*y + (0.00658939)*y + 32.4847; //D1
        if (x>10 && x <20 ) response[1]=   (0.00301674)*y*y + (-0.446544)*y + 27.6374; //D4
        if (x>20 && x <32.8 ) response[2]= (-0.000275694)*y*y + (0.00124251)*y + 18.8999; //D3
        if (x>32.8 && x <42.8 ) response[3]= (-0.00139525)*y*y + (0.104993)*y + 18.1047; //D2
        
            // cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
      //  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    else if (channel==3)// Side Upstream
    {// Assuming an overall size of 42.8 cm with 4 bars of
        double x=-xx/10;
        double y=(yy+346./2)/10.;
        
        double parm[4]={-0.04, -0.05, 1.4, 85.};
           
        
        for(unsigned int s=0; s<4; s++) response[s] =0.;
        response[0] = parm[0]*x*x + parm[1]*y*y + parm[2]*y + parm[3];

        // cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
         // cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    else if (channel==4)// Side Downstream
    {// Assuming an overall size of 42.8 cm with 4 bars of
        double x=xx/10;
        double y=(yy+346./2)/10.;
        
        double parm[4]={-0.04, -0.05, 1.4, 75.};  
        
        
        for(unsigned int s=0; s<4; s++) response[s] =0.;
        response[0] = parm[0]*x*x + parm[1]*y*y + parm[2]*y + parm[3];
        
        //cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
        //cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }

    else if (channel==5)//Right
    {
        double x=-yy/10.;
        double y=(1058/2.-zz)/10.;

        double parm[4][8]={
            {2.34524e+01,	4.28317e-02,	-5.91894e-01,	-5.13309e-03,	-2.47905e-04,	-3.44887e-03,	4.25481e-05,	-1.03817e-05},
            {1.68313e+01, 	5.36853e-02,  	-2.14037e-01,	-4.80535e-03, 	-4.65364e-04,	-1.66572e-02,	4.89028e-05,	-3.33380e-05},
            {2.50310e+01,	-3.10007e-02,	3.57657e-01,	-1.39833e-02,	2.99406e-04,	-3.23669e-02,	1.27237e-04,	-1.13100e-05},
            {1.74834e+01,	1.83925e-01,	5.36737e-01,	7.09769e-04,	-1.64490e-03,	7.48199e-03,	3.43011e-08,	2.11894e-05}
            
        };
        
        for(unsigned int s=0; s<4; s++)
        response[s] = parm[s][7]*x*x*y + parm[s][6]*x*y*y + parm[s][5]*x*x + parm[s][4]*y*y + parm[s][3]*x*y + parm[s][2]*x + parm[s][1]*y + parm[s][0];
        
        //cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
        //cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
       
    }
    else if (channel==6)//Left
    {
        double x=-yy/10.;
        double y=(1058/2.-zz)/10.;
        
        double parm[4][8]={
            {8.12418e+00,	6.61315e-02,	-2.99641e-01,	-9.10408e-04,	-6.79474e-04,	2.00648e-03,	1.24963e-05,	-1.73809e-05},
            {1.19501e+01,	4.76291e-02,	-1.77047e-01,	9.27111e-05,	-4.63061e-04,	-1.40014e-02,	4.39766e-06,	-2.93896e-05},
            {1.68607e+01,	-4.15476e-02,	2.54857e-01,	-6.87363e-03,	3.26876e-04,	-2.65178e-02,	5.62748e-05,	-3.56067e-06},
            {9.73394e+00,	1.56111e-01,	3.27558e-01,	2.45041e-03,	-1.31615e-03,	5.82688e-03,	-1.48528e-05,	2.35177e-05}
            
        };
        
        for(unsigned int s=0; s<4; s++)
        response[s] = parm[s][7]*x*x*y + parm[s][6]*x*y*y + parm[s][5]*x*x + parm[s][4]*y*y + parm[s][3]*x*y + parm[s][2]*x + parm[s][1]*y + parm[s][0];
        
    }

  //  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
    // cout <<  " res[0]: " << response[0] << " channel " << channel<<endl;
    return response;
}

double* veto_HitProcess::IVresponseProposal(int channel, double xx, double yy,double zz, double sx, double sy, double sz)
{
    // Response of the different IV plastic paddles
    // ch
    //
    //
    static double response[4];
    
    for(unsigned int s=0; s<4; s++)response[s] = 0.;
    
    if (channel==1)//top
    {
        double x=xx/10.;
        double y=(sz-zz)/10.;
        
        double parm[4][8]={
            {1.99627e+01,	1.64910e-01,	-5.83528e-01,	-7.34483e-03,	-1.25062e-03, 	4.43805e-03,	5.63766e-05, 	1.40682e-05},
            {1.86162e+01,	4.36475e-02,	-6.78752e-02,	-5.47887e-03,	-1.60512e-04,	-2.33958e-02,	5.55285e-05,	-5.94424e-05},
            {1.85966e+01,	1.96301e-01,	1.34868e-01,	-7.66131e-04,	-1.61720e-03,	-1.91598e-02,	-1.76198e-06,	-4.72970e-05},
            {9.73394e+00,	1.56111e-01,	3.27558e-01,	2.45041e-03,	-1.31615e-03,	5.82688e-03,	-1.48528e-05,	2.35177e-05}
            
        };
        
        for(unsigned int s=0; s<4; s++)
        response[s] = parm[s][7]*x*x*y + parm[s][6]*x*y*y + parm[s][5]*x*x + parm[s][4]*y*y + parm[s][3]*x*y + parm[s][2]*x + parm[s][1]*y + parm[s][0];
        
        //     cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
        //cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    else if (channel==2)//bottom
    {// Assuming an overall size of 42.8 cm with 4 bars of
        double x=-(xx-sx)/10;
        double y=(sz-zz)/10.;
        
        
        for(unsigned int s=0; s<4; s++) response[s] =0.;
        if (x<10)           response[0]= (-0.000303034)*y*y + (0.00658939)*y + 32.4847; //D1
        if (x>10 && x <20 ) response[1]=   (0.00301674)*y*y + (-0.446544)*y + 27.6374; //D4
        if (x>20 && x <32.8 ) response[2]= (-0.000275694)*y*y + (0.00124251)*y + 18.8999; //D3
        if (x>32.8 && x <42.8 ) response[3]= (-0.00139525)*y*y + (0.104993)*y + 18.1047; //D2
        
        // cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
        //  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    else if (channel==6)// Side Upstream
    {// Assuming an overall size of 42.8 cm with 4 bars of
        double x=-xx/10;
        double y=(yy+sy)/10.;
        
        double parm[4]={-0.04, -0.05, 1.4, 85.};
        
        
        for(unsigned int s=0; s<4; s++) response[s] =0.;
        response[0] = parm[0]*x*x + parm[1]*y*y + parm[2]*y + parm[3];
        
        // cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
        // cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    else if (channel==5)// Side Downstream
    {// Assuming an overall size of 42.8 cm with 4 bars of
        double x=xx/10;
        double y=(yy+sy)/10.;
        
        double parm[4]={-0.04, -0.05, 1.4, 75.};
        
        
        for(unsigned int s=0; s<4; s++) response[s] =0.;
        response[0] = parm[0]*x*x + parm[1]*y*y + parm[2]*y + parm[3];
        
        //cout <<  " x: " << x <<  " y: " << y << " yy: " << yy << endl ;
        //cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    
    else if (channel==3)//Right
    {
        double x=-yy/10.;
        double y=(sz-zz)/10.;
 
        double parm[4][8]={
            {2.34524e+01,	4.28317e-02,	-5.91894e-01,	-5.13309e-03,	-2.47905e-04,	-3.44887e-03,	4.25481e-05,	-1.03817e-05},
            {1.68313e+01, 	5.36853e-02,  	-2.14037e-01,	-4.80535e-03, 	-4.65364e-04,	-1.66572e-02,	4.89028e-05,	-3.33380e-05},
            {2.50310e+01,	-3.10007e-02,	3.57657e-01,	-1.39833e-02,	2.99406e-04,	-3.23669e-02,	1.27237e-04,	-1.13100e-05},
            {1.74834e+01,	1.83925e-01,	5.36737e-01,	7.09769e-04,	-1.64490e-03,	7.48199e-03,	3.43011e-08,	2.11894e-05}
            
        };
        
        for(unsigned int s=0; s<4; s++)
        response[s] = parm[s][7]*x*x*y + parm[s][6]*x*y*y + parm[s][5]*x*x + parm[s][4]*y*y + parm[s][3]*x*y + parm[s][2]*x + parm[s][1]*y + parm[s][0];
        
        //cout <<  " x: " << x <<  " y: " << y << " zz: " << zz << endl ;
        //cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
        
    }
    else if (channel==4)//Left
    {
        double x=-yy/10.;
        double y=(sz-zz)/10.;
        
        double parm[4][8]={
            {8.12418e+00,	6.61315e-02,	-2.99641e-01,	-9.10408e-04,	-6.79474e-04,	2.00648e-03,	1.24963e-05,	-1.73809e-05},
            {1.19501e+01,	4.76291e-02,	-1.77047e-01,	9.27111e-05,	-4.63061e-04,	-1.40014e-02,	4.39766e-06,	-2.93896e-05},
            {1.68607e+01,	-4.15476e-02,	2.54857e-01,	-6.87363e-03,	3.26876e-04,	-2.65178e-02,	5.62748e-05,	-3.56067e-06},
            {9.73394e+00,	1.56111e-01,	3.27558e-01,	2.45041e-03,	-1.31615e-03,	5.82688e-03,	-1.48528e-05,	2.35177e-05}
            
        };
        
        for(unsigned int s=0; s<4; s++)
        response[s] = parm[s][7]*x*x*y + parm[s][6]*x*y*y + parm[s][5]*x*x + parm[s][4]*y*y + parm[s][3]*x*y + parm[s][2]*x + parm[s][1]*y + parm[s][0];
        
    }
    
    //  cout <<  " res[0]: " << response[0] <<  " res[1]: " << response[1]<<  " res[2]: " << response[2]<<  " res[3]: " << response[3]  << endl ;
    // cout <<  " res[0]: " << response[0] << " channel " << channel<<endl;
    return response;
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



// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> veto_HitProcess :: electronicNoise()
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
map< int, vector <double> > veto_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double veto_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}




