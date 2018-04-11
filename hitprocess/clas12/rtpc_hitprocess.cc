#include "rtpc_hitprocess.h"
#include <math.h>
#include <random>

// G4 Headers
#include "Randomize.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


map<string, double> rtpc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
    // Establish constants
    PAD_W =2.79;
    PAD_L = 4.0;
    PAD_S = 80.0;
    RTPC_L = 384.0;
    phi_per_pad=PAD_W/PAD_S;
    
    // parameters to change for gas mixture and potential
    // gas mixture = He:CO2 at 80:20
    // potential = 3500V
    a_t=1741.179712;
    b_t=-1.25E+02;
    c_t=388.7449859;
    d_t=-4.33E+01;
    
    a_phi=0.161689123;
    b_phi=0.023505021;
    c_phi=6.00E-06;
    d_phi=2.00E-06;
    
    a_z=0.035972097;
    b_z =-0.000739386;
    
    t_2GEM2 = 296.082;
    sigma_t_2GEM2 = 8.72728;
    t_2GEM3 = 296.131;
    sigma_t_2GEM3 = 6.77807;
    t_2PAD = 399.09;
    sigma_t_2PAD = 7.58056;
    t_2END = t_2GEM2 + t_2GEM3 + t_2PAD;
    sigma_t_gap = sqrt(pow(sigma_t_2GEM2,2)+pow(sigma_t_2GEM3,2)+pow(sigma_t_2PAD,2));
    
    phi_2GEM2 = 0.0492538;
    sigma_phi_2GEM2 = 0.00384579;
    phi_2GEM3 = 0.0470817;
    sigma_phi_2GEM3 = 0.00234478;
    phi_2PAD = 0.0612122;
    sigma_phi_2PAD = 0.00238653;
    phi_2END = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
    sigma_phi_gap = sqrt(pow(sigma_phi_2GEM2,2)+pow(sigma_phi_2GEM3,2)+pow(sigma_phi_2PAD,2));
    
    TPC_TZERO = 0.0;
    
     int chan=0;

    
    static const double PI=3.1415926535;
    
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// true information
	// for example tInfos.eTot is total energy deposited
	trueInfos tInfos(aHit);
	
	// local variable for each step
	vector<G4ThreeVector> Lpos = aHit->GetLPos();

	// take momentum for each step
	vector<G4ThreeVector> Lmom = aHit->GetMoms();

	// energy at each step
	// so tInfos.eTot is the sum of all steps s of Edep[s] 
	vector<double>      Edep = aHit->GetEdep();

	//
	// Get the information x,y,z and Edep at each ionization point
	// 

	double LposX=0.;
	double LposY=0.;
	double LposZ=0.;
	double DiffEdep=0.;
    
	if(tInfos.eTot > 0)

	  {
	    
	    int adc =0;
	    double tdc =0;
          
          
	    for(unsigned int s=0; s<tInfos.nsteps; s++)
	    {
	      LposX = Lpos[s].x();
	      LposY = Lpos[s].y();
	      LposZ = Lpos[s].z();

              if(sqrt(LposX*LposX+LposY*LposY)<30. || sqrt(LposX*LposX+LposY*LposY)>70.) continue;
              if(LposZ<-200. || LposZ>200.) continue;
            
	      DiffEdep = Edep[s];


	      double r0,phi0_rad;
	      //convert (x0,y0,z0) into (r0,phi0,z0)
          r0=(sqrt(LposX*LposX+LposY*LposY))/10.0;  //in cm
          
            phi0_rad=atan2(LposY,LposX); //return (-Pi, + Pi)
	      if( phi0_rad<0.)  phi0_rad+=2.0*PI;
          if( phi0_rad>=2.0*PI )  phi0_rad-=2.0*PI;

          
          // -------------------------------- Addition of Diffusion -----------------------------
          
          // calculate drift time [ns] to first GEM
            double t_drift = a_t*(7.0-r0)+b_t*(7.0-r0)*(7.0-r0);
            
          // determine sigma of drift time [ns]
           double t_diff = sqrt(c_t*(7.0-r0)+d_t*(7.0-r0)*(7.0-r0));

          // calculate drift angle to first GEM at 7 cm [rad]
          double phi_drift = a_phi*(7.0-r0)+b_phi*(7.0-r0)*(7.0-r0);
            
          // determine sigma of drift angle [rad]
          double phi_diff = sqrt(c_phi*(7.0-r0)+d_phi*(7.0-r0)*(7.0-r0));

            // calculate drift in z [mm]
            double z_drift = 0.0;
            
          // determine sigma in z [mm]
          double z_diff = sqrt(a_z*(7.0-r0)+b_z*(7.0-r0)*(7.0-r0));
            
            // find t_s2pad and delta_phi by gaussians
            double t_s2pad = G4RandGauss::shoot(t_drift, t_diff);
            double delta_phi = G4RandGauss::shoot(phi_drift, phi_diff);
            double delta_z = G4RandGauss::shoot(z_drift, z_diff);
            double t_gap = G4RandGauss::shoot(t_2END, sigma_t_gap);
            double phi_gap = G4RandGauss::shoot(phi_2END, sigma_phi_gap);
            
          // ------------------------------------------------------------------------------------
	      
            double phi_rad= phi0_rad+delta_phi+phi_gap;   //phi at pad pcb board
	      if( phi_rad<0. )  phi_rad+=2.0*PI;
	      if( phi_rad>2.0*PI )  phi_rad-=2.0*PI;

            tdc=t_s2pad+t_gap;
            adc=DiffEdep;

	    double z_pos = LposZ+delta_z;
            
            int row = (int) (phi_rad/phi_per_pad);
            float z_shift = row%4;

            if( z_pos+RTPC_L/2.0 < z_shift ){
              continue;
            }
            int col = (int) ((z_pos-z_shift+RTPC_L/2.0)/PAD_L);
            int Num_of_Col = (int) (RTPC_L/PAD_L);
            chan = row*Num_of_Col+col;


          dgtz["CellID"] = (int) chan;
          dgtz["ADC"]    = (int) adc;
          dgtz["Time"]   = tdc;
          dgtz["hitn"]   = (int) hitn;
    
	    } // end step

	  }	      

	return dgtz;
}



vector<identifier>  rtpc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
    static const double PI=3.1415926535;
    double r0 = 0.;
    double phi0_rad=0.;
    
    // Establish constants
    PAD_W =2.79;
    PAD_L = 4.0;
    PAD_S = 80.0;
    RTPC_L = 384.0;
    phi_per_pad=PAD_W/PAD_S;
    
    // parameters to change for gas mixture and potential
    // gas mixture = He:CO2 at 80:20
    // potential = 3500V
    a_t=1741.179712;
    b_t=-1.25E+02;
    c_t=388.7449859;
    d_t=-4.33E+01;
    
    a_phi=0.161689123;
    b_phi=0.023505021;
    c_phi=6.00E-06;
    d_phi=2.00E-06;
    
    a_z=0.035972097;
    b_z =-0.000739386;
    
    phi_2GEM2 = 0.0492538;
    sigma_phi_2GEM2 = 0.00384579;
    phi_2GEM3 = 0.0470817;
    sigma_phi_2GEM3 = 0.00234478;
    phi_2PAD = 0.0612122;
    sigma_phi_2PAD = 0.00238653;
    phi_2END = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
    sigma_phi_gap = sqrt(pow(sigma_phi_2GEM2,2)+pow(sigma_phi_2GEM3,2)+pow(sigma_phi_2PAD,2));

  vector<identifier> yid = id;

  G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();
  G4ThreeVector Lxyz = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()
    ->GetTopTransform().TransformPoint(xyz);///< Local Coordinates of interaction

    int chan;
  double LposX = 0;
  double LposY = 0;
  double LposZ = 0;
  LposX = Lxyz.x();
  LposY = Lxyz.y();
  LposZ = Lxyz.z();
    
  r0=(sqrt(LposX*LposX+LposY*LposY))/10.0;  //in mm
  
  phi0_rad=atan2(LposY,LposX); //return (-Pi, + Pi)
  if( phi0_rad<0.)  phi0_rad+=2.0*PI;
  if( phi0_rad>=2.0*PI)  phi0_rad-=2.0*PI;
  
    // -------------------------------- Addition of Diffusion -----------------------------
    
    // calculate drift angle to first GEM at 7 cm [rad]
    double phi_drift = a_phi*(7.0-r0)+b_phi*(7.0-r0)*(7.0-r0);
    
    // determine sigma of drift angle [rad]
    double phi_diff = sqrt(c_phi*(7.0-r0)+d_phi*(7.0-r0)*(7.0-r0));
    
    // calculate drift in z [mm]
    double z_drift = 0.0;
    
    // determine sigma in z [mm]
    double z_diff = sqrt(a_z*(7.0-r0)+b_z*(7.0-r0)*(7.0-r0));
    
    // find t_s2pad and delta_phi by gaussians
    double delta_phi = G4RandGauss::shoot(phi_drift, phi_diff);
    double delta_z = G4RandGauss::shoot(z_drift, z_diff);
    double phi_gap = G4RandGauss::shoot(phi_2END, sigma_phi_gap);

    
    // ------------------------------------------------------------------------------------
    
    double phi_rad= phi0_rad+delta_phi+phi_gap;   //phi at pad pcb board
	      if( phi_rad<0. )  phi_rad+=2.0*PI;
	      if( phi_rad>2.0*PI )  phi_rad-=2.0*PI;

	double z_pos = LposZ+delta_z;
    
    int row = int(phi_rad/phi_per_pad);
    float z_shift = row%4;
    int col = (int) ((z_pos-z_shift+RTPC_L/2.0)/PAD_L);
    int Num_of_Col = (int) (ceil(RTPC_L/PAD_L));
    chan = row*Num_of_Col+col;

    yid[0].id=chan;
  return yid;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> rtpc_HitProcess :: electronicNoise()
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


map< string, vector <int> >  rtpc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > rtpc_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double rtpc_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}
