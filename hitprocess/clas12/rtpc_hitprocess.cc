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
	static const double PI=3.1415926535;
	rejectHitConditions = false;
	writeHit = true;
	
	// Establish constants
	// NEW CONSTANTS FOR UPDATED RECONSTRUCTION
	PAD_W =2.79;
	PAD_L = 4.0;
	PAD_S = 79.0;//=80.0 old version after comment
	RTPC_L = 384.0;
	phi_per_pad=(2.0*PI)/180;
	
	// Drift time from first GEM to readout pad
	t_2GEM2 = 0.0; // 169.183;
	sigma_t_2GEM2 = 8.72728;
	t_2GEM3 = 0.0; //222.415;
	sigma_t_2GEM3 = 5.62223;
	t_2PAD = 0.0; //414.459;
	sigma_t_2PAD = 7.58056;
	
	phi_2GEM2 = 0.0; // 0.0416925;
	sigma_phi_2GEM2 = 0.00384579;
	phi_2GEM3 = 0.0; // 0.0416574;
	sigma_phi_2GEM3 = 0.00160235;
	phi_2PAD = 0.0; // 0.057566;
	sigma_phi_2PAD = 0.00238653;
	
	// parameters to change for gas mixture and potential
	// gas mixture = He:CO2 at 85.5:14.5
	// potential = 4450V
	
	// --------- Updated 25 June 2020 B-field parameters ----------- //
	double a_t0 = 1470.0; // -0.607056;
	double a_t1 = 0.009; // -0.00138137;
	//double a_t2 = -9.0e05; // 1.6494e-05; //Gave crazy big -ve values to tdc or time 
	double a_t2 = -9.0e-05; // 1.6494e-05;
	
	double b_t0 = 3290.0; // 3183.62;
	double b_t1 = -0.0138; // -0.0269385;
	double b_t2 = -0.000332; // -0.00037645;
	
	double c_t0 = 0.0; // -0.0200932;
	double c_t1 = 0.0; // 0.00319441;
	double c_t2 = 0.0; // -8.75299e-06;
	
	double a_phi0 = -0.113; // 0;
	double a_phi1 = -1.05e-6; // 0;
	double a_phi2 = 1.58e-8; // 0;
	
	double b_phi0 = -0.9697; // 1.22249;
	double b_phi1 = 9.48e-6; // -2.3708e-05;
	double b_phi2 = 9.9e-8; // -1.11055e-07;
	
	double c_phi0 = 0.0; // -1.13197;
	double c_phi1 = 0.0; // 0.000118575;
	double c_phi2 = 0.0; // 2.50094e-08;
	
	
	// Diffusion parameters
	diff_at = 388.7449859;
	diff_bt = -4.33E+01;
	
	diff_aphi = 6.00E-06;
	diff_bphi = 2.00E-06;
	
	a_z=0.035972097;
	b_z =-0.000739386;
	
	t_2END = t_2GEM2 + t_2GEM3 + t_2PAD;
	sigma_t_gap = sqrt(pow(sigma_t_2GEM2,2)+pow(sigma_t_2GEM3,2)+pow(sigma_t_2PAD,2));
	phi_2END = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
	sigma_phi_gap = sqrt(pow(sigma_phi_2GEM2,2)+pow(sigma_phi_2GEM3,2)+pow(sigma_phi_2PAD,2));
	// ----------------------------------------- //
	
	// ------------ B=0 T Parameters ------------//
	/*double a_t1 = 0;
	 double a_t2 = 0;
	 double a_t3 = 0;
	 double a_t4 = 0;
	 double a_t5 = 6.96387e+02;
	 
	 
	 double b_t1 = 0;
	 double b_t2 = 0;
	 double b_t3 = 0;
	 double b_t4 = 0;
	 double b_t5 = -4.73759e+01;
	 
	 double a_phi1 = 0;
	 double a_phi2 = 0;
	 double a_phi3 = 0;
	 double a_phi4 = 0;
	 double a_phi5 = 0;
	 
	 double b_phi1 = 0;
	 double b_phi2 = 0;
	 double b_phi3 = 0;
	 double b_phi4 = 0;
	 double b_phi5 = 0;
	 
	 // Diffusion parameters
	 c_t=0.0;
	 d_t=0.0;
	 
	 c_phi=0.0;
	 d_phi=0.0;
	 
	 a_z=0.0;
	 b_z =0.0;
	 
	 t_2END = 500.0;
	 sigma_t_gap = 0;
	 phi_2END = 0.0;
	 sigma_phi_gap = 0.0;*/
	// ------------------------------------------//
	
	
	TPC_TZERO = 0.0;
	
	// not used anymore?
	//int chan=0;
	
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
	
	
	// -------------------------- TIME SHIFT for non-primary tracks ---------------------------
	//M. U.: map.find(Key) tries to find an entry with that key, if it doesn't find it, it will return map.end()
	//     So, to ensure that hits of a given track gets the same timeShift (but random among different tracks)
	//     we must check if the find returns map.end(). If it does, we should create a new entry, else we should 
	//     skip from overwriting/reassigning in order to keep the same old value for all hits of a given track.
	//if(aHit->GetTId() != timeShift_map.cend()->first){
	if(timeShift_map.find(aHit->GetTId()) == timeShift_map.end()){
		if(aHit->GetTId()< 3) timeShift_map[aHit->GetTId()] = 0.0;
		else timeShift_map[aHit->GetTId()] = G4RandFlat::shoot(-8000.,8000.);
	}
	
	//cout<<"kp: aHit->GetTId(): "<<aHit->GetTId()<<endl;
	vector<int> tids = aHit->GetTIds();
	vector<int> mtids = aHit->GetmTrackIds();  
	//raws["otid"]    = (double) aHit->GetoTrackId(); 
	vector<int>    otids = aHit->GetoTrackIds();
	/*
	cout<<"kp: tInfos.nsteps = "<<tInfos.nsteps<<" aHit->GetTIds().size() = "<<tids.size()
	    <<" aHit->GetmTrackIds().size() = "<<mtids.size()<<" aHit->GetoTrackIDs().size() = "<<otids.size()<<endl; 
	cout<<"tids[step]   mtids[step]   otids[step]"<<endl;
	for (unsigned int s=0; s<tInfos.nsteps; s++) cout<< tids[s] << "   "<< mtids[s] << "   "<< otids[s]<<endl;
	*/
	
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
			
			//          if(sqrt(LposX*LposX+LposY*LposY)<30. || sqrt(LposX*LposX+LposY*LposY)>70.) continue;
			//          if(LposZ<-(RTPC_L/2.0) || LposZ>(RTPC_L/2.0)) continue;
			
			DiffEdep = Edep[s];
			
			z_cm = LposZ/10.0;
			
			// all in cm
			a_t= a_t0 + a_t1*(pow(z_cm,2)) + a_t2*(pow(z_cm,4));
			b_t= b_t0 + b_t1*(pow(z_cm,2)) + b_t2*(pow(z_cm,4));
			c_t= c_t0 + c_t1*(pow(z_cm,2)) + c_t2*(pow(z_cm,4));
			     
			a_phi=a_phi2*(pow(z_cm,4)) + a_phi1*(pow(z_cm,2)) + a_phi0;
			b_phi=b_phi2*(pow(z_cm,4)) + b_phi1*(pow(z_cm,2)) + b_phi0;
			c_phi=c_phi2*(pow(z_cm,4)) + c_phi1*(pow(z_cm,2)) + c_phi0;
			
			double r0,phi0_rad;
			//convert (x0,y0,z0) into (r0,phi0,z0)
			r0=(sqrt(LposX*LposX+LposY*LposY))/10.0;  //in cm
			
			phi0_rad=atan2(LposY,LposX); //return (-Pi, + Pi)
			if( phi0_rad<0.)  phi0_rad+=2.0*PI;
			if( phi0_rad>=2.0*PI )  phi0_rad-=2.0*PI;
			
			
			// -------------------------------- Addition of Diffusion -----------------------------
			
			// calculate drift time [ns] to first GEM
			double t_drift = a_t + b_t*(7.0*7.0-r0*r0)/40.0;
			
			// determine sigma of drift time [ns]
			double t_diff = sqrt(diff_at*(7.0-r0)+diff_bt*(7.0-r0)*(7.0-r0));
			
			// calculate drift angle to first GEM at 7 cm [rad]
			double phi_drift = a_phi + b_phi*log(7.0/r0) + c_phi*((1.0/(r0*r0)) - (1.0/(49.0)));
			//double phi_drift = a_phi*(7.0-r0)+b_phi*(7.0-r0)*(7.0-r0);
			
			// determine sigma of drift angle [rad]
			double phi_diff = sqrt(diff_aphi*(7.0-r0)+diff_bphi*(7.0-r0)*(7.0-r0));
			
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
			if( phi_rad<0.0 )  phi_rad+=2.0*PI;
			if( phi_rad>=2.0*PI )  phi_rad-=2.0*PI;
			
			// time shift
			//shift_t = timeShift_map.find(aHit->GetTId())->second;
			shift_t = timeShift_map.find(otids[s])->second; //To ensure secondaries have the same shift_t as primaries or acencestors
			
			// NO time shift
			//shift_t = 0.0;
			
			tdc=t_s2pad+t_gap+shift_t; //kp:
			adc=((int)100000.0*DiffEdep); //cludge for tiny ADC numbers
			//adc = 0.5 + adc; //kp: For debugging purpose to avoid lots of zeros (due to int typecasting below)
			//adc = 100000*DiffEdep;//kp: For same reason as bove 
			
			//double z_pos = (z_cm+delta_z)*10.0;// must be mm
			double z_pos = z_cm*10.0 + delta_z;// must be mm
			int col = -999;
			int row = -999;
			
			row = ceil(phi_rad/phi_per_pad);
			float z_shift = row%4;
			
			float col_min = -RTPC_L/2.0+z_shift;
			float col_max = RTPC_L/2.0+z_shift;
			
			if( z_pos < col_min || z_pos > col_max ) {row = -999; col = -999;}
			else {
				
				row = ceil(phi_rad/phi_per_pad);
				col = ceil((z_pos+RTPC_L/2.0-z_shift)/PAD_L);
			}
			
			//int Num_of_Col = (int) (RTPC_L/PAD_L);
			//chan = row*Num_of_Col+col;
			
			//Mohammad: kill all the hits in the last three rows that correspond to the seam effect
			if(row>177){row = -999; col = -999;}
			
			dgtz["Sector"]    = 1;
			dgtz["Layer"]     = col;
			dgtz["Component"] = row;
			dgtz["Order"]     = 0;
			dgtz["Time"]      = tdc;
			dgtz["ADC"]       = (int) adc;
			dgtz["Ped"]       = 0;
			//dgtz["TimeShift"] = shift_t;
			dgtz["hitn"]      = (int) hitn;
			
		} // end step
		
	}
	
	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}



vector<identifier>  rtpc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	static const double PI=3.1415926535;
	
	// Establish constants
	PAD_W = 2.758; //2.79 !!
	PAD_L = 4.0;
	PAD_S = 79.0;
	RTPC_L = 384.0;
	phi_per_pad=PAD_W/PAD_S;
	
	// Drift time from first GEM to readout pad
	phi_2GEM2 = 0.0; // 0.0416925;
	sigma_phi_2GEM2 = 0.00384579;
	phi_2GEM3 = 0.0; // 0.0416574;
	sigma_phi_2GEM3 = 0.00160235;
	phi_2PAD = 0.0; // 0.057566;
	sigma_phi_2PAD = 0.00238653;
	
	// parameters to change for gas mixture and potential
	// gas mixture = He:CO2 at 85.5:14.5
	// potential = 4450V
	
	// --------- Updated 14 Dec 2020 B-field parameters
	double a_phi0 = -0.113; // 0;
	double a_phi1 = -1.05e-6; // 0;
	double a_phi2 = 1.58e-8; // 0;
	double b_phi0 = -0.9697; // 1.22249;
	double b_phi1 = 9.48e-6; // -2.3708e-05;
	double b_phi2 = 9.9e-8; // -1.11055e-07;
	double c_phi0 = 0.0; // -1.13197;
	double c_phi1 = 0.0; // 0.000118575;
	double c_phi2 = 0.0; // 2.50094e-08;

	// Diffusion parameters
	diff_aphi = 6.00E-06;
	diff_bphi = 2.00E-06;
	
	
	a_z = 0.035972097;
	b_z = -0.000739386;
	
	phi_2END = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
	sigma_phi_gap = sqrt(pow(sigma_phi_2GEM2,2)+pow(sigma_phi_2GEM3,2)+pow(sigma_phi_2PAD,2));
	
	vector<identifier> yid = id;
	
	G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector Lxyz = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()
	->GetTopTransform().TransformPoint(xyz);///< Local Coordinates of interaction
	
	int chan = 0;
	double LposX = 0;
	double LposY = 0;
	double LposZ = 0;
	LposX = Lxyz.x();
	LposY = Lxyz.y();
	LposZ = Lxyz.z();
	
	z_cm = LposZ/10.0;
	
	// all in cm
	a_phi=a_phi2*(pow(z_cm,4)) + a_phi1*(pow(z_cm,2)) + a_phi0;
	b_phi=b_phi2*(pow(z_cm,4)) + b_phi1*(pow(z_cm,2)) + b_phi0;
	c_phi=c_phi2*(pow(z_cm,4)) + c_phi1*(pow(z_cm,2)) + c_phi0;
	
	double r0,phi0_rad;
	//convert (x0,y0,z0) into (r0,phi0,z0)
	r0=(sqrt(LposX*LposX+LposY*LposY))/10.0;  //in cm
	
	
	phi0_rad=atan2(LposY,LposX); //return (-Pi, + Pi)
	if( phi0_rad<0.)  phi0_rad+=2.0*PI;
	if( phi0_rad>=2.0*PI)  phi0_rad-=2.0*PI;
	
	// -------------------------------- Addition of Diffusion -----------------------------
	
	// calculate drift angle to first GEM at 7 cm [rad]
	double phi_drift = a_phi + b_phi*log(7.0/r0) + c_phi*((1.0/(r0*r0)) - (1.0/(49.0)));
	//double phi_drift = a_phi*(7.0-r0)+b_phi*(7.0-r0)*(7.0-r0);
	
	// determine sigma of drift angle [rad]
	double phi_diff = sqrt(diff_aphi*(7.0-r0)+diff_bphi*(7.0-r0)*(7.0-r0));
	
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
	if( phi_rad<0.0 )  phi_rad+=2.0*PI;
	if( phi_rad>=2.0*PI )  phi_rad-=2.0*PI;
	
	double z_pos = (z_cm+delta_z)*10.0; //mm
	int col = -999;
	int row = -999;
	
	row = ceil(phi_rad/phi_per_pad);
	float z_shift = row%4;
	
	float col_min = -RTPC_L/2.0+z_shift;
	float col_max = RTPC_L/2.0+z_shift;
	
	if( z_pos < col_min || z_pos > col_max ) {row = -999; col = -999;}
	else {
	
	  row = ceil(phi_rad/phi_per_pad);
	  col = ceil((z_pos+RTPC_L/2.0-z_shift)/PAD_L);
	}
	
	int Num_of_Col = (int) (RTPC_L/PAD_L);
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










