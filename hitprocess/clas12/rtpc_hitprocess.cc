#include "rtpc_hitprocess.h"
#include <math.h>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;


float PAD_W =2.79;
float PAD_L = 4.0;
float PAD_S = 80.0;
float RTPC_L = 400.0;
float phi_per_pad=PAD_W/PAD_S;
float a=-0.177, b=0.0392,c=10.977;
float d=0.027,e=-0.517,f=2.4615;

float TPC_TZERO = 0.0;

static const double PI=3.1415926535;

map<string, double> rtpc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
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
	    int chan=0;
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

          // determine drift time
	      double t_s2pad =((a*r0*r0)+(b*r0)+c)*1000.;  // calculate drift time in ns

          double delta_phi=(d*r0*r0)+(e*r0)+f;   // calculate dphi to edge of RTPC
	      
            double phi_rad= phi0_rad+delta_phi;   //phi at pad pcb board
	      if( phi_rad<0. )  phi_rad+=2.0*PI;
	      if( phi_rad>2.0*PI )  phi_rad-=2.0*PI;

	      tdc=t_s2pad;
	      adc=DiffEdep;
            
            int row = (int) (phi_rad/phi_per_pad);
            float z_shift = row%4;

            if( LposZ+RTPC_L/2.0 < z_shift ){ 
              continue;
            }

//cout << s << "  " << Edep[s] << "  " << LposX << "  " << LposY << "  " << LposZ << endl;

            z_shift = 0.;
            int col = (int) ((LposZ-z_shift+RTPC_L/2.0)/PAD_L);
            int Num_of_Col = (int) (RTPC_L/PAD_L);
            chan = row*Num_of_Col+col;

//if(LposZ<PAD_L || LpoZ>RTPC_L/2.0) cout << LposZ << "  " << col << "  " << chan << endl;
            
	      //chan = identity[0].id;

          dgtz["phiRad"]   = phi_rad;
	  dgtz["CellID"] = chan;
	  dgtz["ADC"]    = adc;
	  dgtz["TDC"]    = tdc;
          dgtz["PosX"]    = LposX;
          dgtz["PosY"]    = LposY;
          dgtz["PosZ"]    = LposZ;
          dgtz["EDep"]   = DiffEdep;
	  dgtz["step"]   = s;
          dgtz["hitn"]   = hitn;
    
	    } // end step

	  }	      

	return dgtz;
}



vector<identifier>  rtpc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
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
  
  double r0,phi0_rad;
    
  r0=(sqrt(LposX*LposX+LposY*LposY))/10.0;  //in mm
  
  phi0_rad=atan2(LposY,LposX); //return (-Pi, + Pi)
  if( phi0_rad<0.)  phi0_rad+=2.0*PI;
  if( phi0_rad>=2.0*PI)  phi0_rad-=2.0*PI;
  
    double delta_phi=(d*r0*r0)+(e*r0)+f;   // calculate dphi to edge of RTPC in rad    
    double phi_rad= phi0_rad+delta_phi;   //phi at pad pcb board
    if( phi_rad<0. )  phi_rad+=2.0*PI;
    if( phi_rad>2.0*PI )  phi_rad-=2.0*PI;
    
    int row = int(phi_rad/phi_per_pad);
//    float z_shift = row%4;
    float z_shift = 0.;
    int col = (int) ((LposZ-z_shift+RTPC_L/2.0)/PAD_L);
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














