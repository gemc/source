// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "cnd_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static cndConstants initializeCNDConstants(int runno)
{
	// all these constants should be read from CCDB
	cndConstants cndc;

        cout<<"Entering initializeCNDConstants"<<endl;
	if(runno == -1) return cndc;
	// database
	cndc.runNo = runno;
	cndc.date       = "2016-11-27";
	if(getenv ("CCDB_CONNECTION") != NULL)
		cndc.connection = (string) getenv("CCDB_CONNECTION");
	else
		cndc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	cndc.variation  = "default";

	int isec,ilay,istr;
	
	vector<vector<double> > data;

	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(cndc.connection));
        cout<<"Connecting to "<<cndc.connection<<"/calibration/cnd"<<endl;

	cout<<"CND:Getting status"<<endl;
	sprintf(cndc.database,"/calibration/cnd/status:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  //	  cndc.status[isec-1][ilay-1][0].push_back(data[row][3]);
	  cndc.status[isec-1][ilay-1][istr-1]=data[row][3];
	}
	cout<<"CND:Getting attenuation"<<endl;
	sprintf(cndc.database,"/calibration/cnd/att_length:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  //	  cndc.att_length[isec-1][ilay-1][0].push_back(data[row][3]);
	  cndc.att_length[isec-1][ilay-1][istr-1]=data[row][3];
	}
        cout<<"CND:Getting effective_velocity"<<endl;
	sprintf(cndc.database,"/calibration/cnd/veff:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  //	  cndc.veff[isec-1][ilay-1][0].push_back(data[row][3]);
	  cndc.veff[isec-1][ilay-1][istr-1]=data[row][3];
	}
        cout<<"CND:Getting energy calibration"<<endl;
	sprintf(cndc.database,"/calibration/cnd/ecal:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  cndc.ecalD[isec-1][ilay-1][istr-1]=data[row][3];
	  cndc.ecalN[isec-1][ilay-1][istr-1]=data[row][5];
	}
        cout<<"CND:Getting u-turn eloss"<<endl;
	sprintf(cndc.database,"/calibration/cnd/uturn_eloss:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  cndc.uturn_e[isec-1][ilay-1][0]=data[row][3];
	}
        cout<<"CND:Getting u-turn delay"<<endl;
	sprintf(cndc.database,"/calibration/cnd/uturn_tloss:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr = data[row][2];
	  cndc.uturn_t[isec-1][ilay-1][0]=data[row][3];
	}
        cout<<"CND:Getting time offset LR"<<endl;
	sprintf(cndc.database,"/calibration/cnd/time_offsets_LR:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  cndc.time_offset_LR[isec-1][ilay-1][0]=data[row][3];
	}
        cout<<"CND:Getting time offset layer"<<endl;
	sprintf(cndc.database,"/calibration/cnd/time_offsets_layer:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);	
	for(unsigned row = 0; row < data.size(); row++)
	{
	  isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	  cndc.time_offset_layer[isec-1][ilay-1][0]=data[row][3];
	}

	return cndc;
}

map<string, double> cnd_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	string hd_msg = " > cnd hit process";
	
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);
	
	int sector  = identity[0].id;
	int layer  = identity[1].id;
	int paddle = identity[2].id;
	int paddleN;

	// Energy propagation:
	
	double Lg[3];                                 // lengths along light-guides to PMT (from design drawings)
	Lg[0] = 139.56*cm;                            // inner layer
	Lg[1] = 138.89*cm;
	Lg[2] = 136.88*cm;                            // outer layer
	
	//double att_length = 1.5*m;                    // light attenuation length in scintillator
	double att_length_lg = 950*cm;                 // light attenuation length in long light-guides
	
	double sensor_surface = pow(2.5*cm,2)*pi;     // X-sectional area of PMT, assume radius of 2.5 cm.
	double paddle_xsec = 0.;                       // cross-sectional area of the paddle
	if (layer == 1) paddle_xsec = 22.8*cm*cm;     // from the geometry files
	else if (layer == 2) paddle_xsec = 26.4*cm*cm;
	else if (layer == 3) paddle_xsec = 27.7*cm*cm;
	double light_coll;                            // ratio of photo_sensor area over paddle section, times optical coupling ~ light collection efficiency
	if (sensor_surface < paddle_xsec) light_coll = 0.7 * sensor_surface / paddle_xsec;
	else light_coll = 0.7;                        // to make sure sensor_surface / paddle_xsec doesn't go over 1.
	
	double dEdxMIP = 1.956;  
	double thickness = 3;
	
	// double uturn[3];                               // fraction of energy which makes it through the u-turn light-guides (based on cosmic tests)
	// uturn[0] = 0.65;                               // inner layer
	// uturn[1] = 0.6;
	// uturn[2] = 0.5;                                // outer layer
	
	// double light_yield = 10000;               // number of optical photons produced in the scintillator per MeV of deposited energy
	// double sensor_qe = 0.2;                       // photo sensor quantum efficiency
	// double sensor_gain = 0.24;                    // gain of the photo sensor in pC/(#p.e.); it defines the conversion from photoelectrons to charge:
						      // for a pmt gain of 1.5*10^6, this factor is equal to 1.5*10^6*1.6*10^-19 C = 0.24 pC/(#p.e.)

	double pmtPEYld = 500;

	// double signal_split = 0.5;                    // signal is split into two, going to QDC and discriminators.
	// double adc_conv = 10.;                        // conversion factor from pC to ADC (typical sensitivy of CAEN VME QDC is of 0.1 pC/ch)
	// double adc_ped = 3.;                          // ADC Pedestal
	
	
	// Time of signal:
	
	//	double veff = 16*cm/ns;                       // light velocity in scintillator

	//	double sigmaTD = 0.14*ns/sqrt(MeV);           // time smearing factor (estimated from tests at Orsay), same paddle as hit (in ns/sqrt(MeV)).
	//double sigmaTN = 0.14*ns/sqrt(MeV);           // time smearing factor, neighbouring paddle to the hit one, units as above.
	double sigmaTD = 0.14;
	double sigmaTN = 0.14;

	double tdc_conv = 40;                      // TDC conversion factor (1/0.025ns), channels per per ns.
		
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
	double e_dir_old = 0.;     // attenuated energy as it arrives at the two PMTs (one coupled to the hit paddle, one to its neighbour)
	double e_neigh_old = 0.;
	
	// Variables for each hit:
	
	//  int flag_counted = 0;   // Flag to specify that the time of hit has already been determined for this hit (1: yes, 0: no)
	
	double et_D = 0.;      // variables to hold total energy*time values for each step, for the two PMTs
	double et_N = 0.;
	double et_D_old = 0.;      // variables to hold total energy*time values for each step, for the two PMTs
	double et_N_old = 0.;

	double etotD = 0.;  // total energy of hit propagated to the PMT connected to the hit paddle
	double etotN = 0.;  // total energy of hit propagated to downstream end of the hit paddle, round u-turn, along neighbouring paddle and light-guide and into PMT
	double etotD_old = 0.;  // total energy of hit propagated to the PMT connected to the hit paddle
	double etotN_old = 0.;  // total energy of hit propagated to downstream end of the hit paddle, round u-turn, along neighbouring paddle and light-guide and into PMT
	double timeD = 0.;  // hit times measured at the upstream edges of the two paddles
	double timeN = 0.;
	//double timeD_old = 0.;  // hit times measured at the upstream edges of the two paddles
	//double timeN_old = 0.;
	
	int ADCD = 0;
	int ADCN = 0;
	int TDCD = 16384;    // max value of the ADC readout
	int TDCN = 16384;
	
    int paddle_N = 3 - paddle;
	double attlength_D = cndc.att_length[sector-1][layer-1][paddle-1];
    double attlength_N = cndc.att_length[sector-1][layer-1][paddle_N-1];
    double v_eff_D = cndc.veff[sector-1][layer-1][paddle-1];
    double v_eff_N = cndc.veff[sector-1][layer-1][paddle_N-1];

	double t_u = cndc.uturn_t[sector-1][layer-1][0];
	double uturn = cndc.uturn_e[sector-1][layer-1][0];
	double t_offset_LR = cndc.time_offset_LR[sector-1][layer-1][0];
	double t_offset_layer = cndc.time_offset_layer[sector-1][layer-1][0];
	double adcd_mip = cndc.ecalD[sector-1][layer-1][paddle-1];
	// double adcn_mip = cndc.ecalN[sector-1][layer-1][paddle-1];

	double dUp = length + tInfos.z;
	double dDn = length - tInfos.z;
	double attUp  = exp(-dUp/cm/attlength_D);
	double attDn  = uturn*exp(-dDn/cm/attlength_D)*exp(-2*length/cm/attlength_N);
	double gain = sqrt(attUp*attDn);

        //if(layer==1)
	//{
	//  t_offset_layer=1.2;
	//  t_offset_LR=-0.7;
	//}
        //if(layer==2)
	//{
	//  t_offset_layer=-1.2;
	//  t_offset_LR=+0.8;
	//}
        //if(layer==3)t_offset_LR=1.2;

	if(tInfos.eTot>0)
	{
		for(unsigned int s=0; s<nsteps; s++)
		  {

		    // Distances travelled through the paddles to the upstream edges (of the hit paddle and of its coupled neighbour):
		    dDir   = (length + Lpos[s].z());
		    dNeigh = (length - Lpos[s].z());
		    
		    // apply Birks effect
		    Edep_B = BirksAttenuation(Edep[s],dx[s],charge[s],birks_constant);
		    // Calculate attenuated energy which will reach both PMTs:
		    e_dir   = (Edep_B/2.) * exp(-dDir/cm/attlength_D);
		    e_neigh = (Edep_B/2.) * exp(-dNeigh/cm/attlength_D) * uturn * exp(-2*length/cm/attlength_N);
		    e_dir_old   = (Edep_B/2.) * exp(-dDir/cm/attlength_D-Lg[layer-1]/att_length_lg) * light_coll;
		    e_neigh_old = (Edep_B/2.) * exp(-dNeigh/cm/attlength_D-Lg[layer-1]/att_length_lg) * uturn * light_coll;
		    // Integrate energy over entire hit. These values are the output and will be digitised:
		    etotD = etotD + e_dir;
		    etotN = etotN + e_neigh;
		    etotD_old = etotD_old + e_dir_old;
		    etotN_old = etotN_old + e_neigh_old;

		    //cout << "step: " << s << " etotD, etotN: " << etotD << ", " << etotN << endl;
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

		    // cout<<"times[s] (ns) "<<times[s]/ns<<endl;

		    et_D = et_D + ((times[s] + dDir/cm/v_eff_D) * e_dir);
		    et_N = et_N + ((times[s] + dNeigh/cm/v_eff_D + 2*length/cm/v_eff_N+ t_u) * e_neigh);
		    et_D_old = et_D_old + ((times[s] + dDir/cm/v_eff_D) * e_dir_old);
		    et_N_old = et_N_old + ((times[s] + dNeigh/cm/v_eff_D + t_u) * e_neigh_old);
		    
		  }   // close loop over steps s
		
		
		
		/**** The following calculates the time based on energy-weighted average of all step times ****/

		timeD = et_D / etotD;      // sum(energy*time) /  sum(energy)
		timeN = et_N / etotN;
		//timeD_old = et_D_old / etotD_old;      // sum(energy*time) /  sum(energy)
		//timeN_old = et_N_old / etotN_old;

		timeD=timeD+t_offset_layer;
		timeN=timeN+t_offset_layer;
		
		if(paddle==2)timeD=timeD+t_offset_LR;
		if(paddle==1)timeN=timeN+t_offset_LR;

		/******** end timing determination ***********/
		// cout<< " Length (cm) "<<length/cm<<endl;
		// cout << "Reconstructed time (ns): " << ((timeD + timeN)/2. - 2.*length/v_eff/cm - t_u/2.)/ns << endl;
		// cout << "Reconstructed z (cm, wrt paddle center): " << ((length + t_u*v_eff/2. - (timeN - timeD)*v_eff/2.)/10.)/cm << endl;
		//cout << "etotD, etotN (in MeV): " << etotD/MeV << ", " << etotN/MeV << endl;
		// cout << "timeD, timeN (in ns): " << timeD/ns << ", " << timeN/ns << endl;
		
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
		  TDCD = (int) ((G4RandGauss::shoot(timeD,sigmaTD/sqrt(etotD)))* tdc_conv);
			// int ADCD_old = (int) (G4Poisson(etotD_old*light_yield*sensor_qe)*signal_split*sensor_gain*adc_conv + adc_ped);
			double npheD = G4Poisson(etotD*pmtPEYld);
			double eneD = npheD/pmtPEYld;
			ADCD = eneD*(adcd_mip/(dEdxMIP*thickness)/gain);
		}
		if (etotN > 0.)
		{
		  TDCN = (int) ((G4RandGauss::shoot(timeN,sigmaTN/sqrt(etotN)))* tdc_conv);
			//TDCN = (int) ((timeN_old + G4RandGauss::shoot(0.,sigmaTN/sqrt(etotN_old*MeV))) * tdc_conv);
			// int ADCN_old = (int) (G4Poisson(etotN_old*light_yield*sensor_qe)*signal_split*sensor_gain*adc_conv + adc_ped);
			double npheN = G4Poisson(etotN*pmtPEYld);
			double eneN = npheN/pmtPEYld;
			ADCN = eneN*(adcd_mip/(dEdxMIP*thickness)/gain);
		}
		
		//if(TDCD < 0) TDCD = 0;
		if (TDCD > 16384) TDCD = 16384;
		//if(TDCN < 0) TDCN = 0;
		if(TDCN > 16384) TDCN = 16384;
		
		if(ADCD < 0) ADCD = 0;
		if(ADCN < 0) ADCN = 0;
		
	}  // closes tInfos.eTot>0
	
	
	if(verbosity>4)
	{
		cout <<  hd_msg << " layer: " << layer    << ", paddle: " << paddle  << " x=" << tInfos.x/cm << "cm, y=" << tInfos.y/cm << "cm, z=" << tInfos.z/cm << "cm" << endl;
		cout <<  hd_msg << " Etot=" << tInfos.eTot/MeV     << "MeV, average time=" << tInfos.time  << "ns"  << endl;
		cout <<  hd_msg << " TDCD= " << TDCD     << ", TDCN= " << TDCN    << ", ADCD= " << ADCD << ", ADCN= " << ADCN << endl;
	}
	if(paddle==1)paddleN=2;
	if(paddle==2)paddleN=1;
	// Status flags
	switch (cndc.status[sector-1][layer-1][paddle-1])
	{
		case 0:
			break;
		case 1:
			ADCD = 0;
			break;
		case 2:
			TDCD = 0;
			break;
		case 3:
			ADCD = TDCD = 0;
			break;

		case 5:
			break;
			
		default:
			cout << " > Unknown CND status: " << cndc.status[sector-1][layer-1][paddle-1] << " for sector " << sector << ",  layer " << layer << ", paddle " << paddle << endl;
	}
	switch (cndc.status[sector-1][layer-1][paddleN-1])
	{
		case 0:
			break;
		case 1:
			ADCN = 0;
			break;
		case 2:
			TDCN = 0;
			break;
		case 3:
			ADCN = TDCN = 0;
			break;

		case 5:
			break;
			
		default:
			cout << " > Unknown CND status: " << cndc.status[sector-1][layer-1][paddleN-1] << " for sector " << sector << ",  layer " << layer << ", paddleN " << paddleN << endl;
	}

	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["layer"]  = layer;
	dgtz["paddle"] = paddle;
	dgtz["ADCD"]   = (int) ADCD;
	dgtz["ADCN"]   = (int) ADCN;
	dgtz["TDCD"]   = (int) TDCD;
	dgtz["TDCN"]   = (int) TDCN;
	
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

// - charge: returns charge/time digitized information / step
map< int, vector <double> > cnd_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double cnd_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}


void cnd_HitProcess::initWithRunNumber(int runno)
{
	if(cndc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		cndc = initializeCNDConstants(runno);
		cndc.runNo = runno;
	}
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> cnd_HitProcess :: electronicNoise()
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

// this static function will be loaded first thing by the executable
cndConstants cnd_HitProcess::cndc = initializeCNDConstants(-1);





