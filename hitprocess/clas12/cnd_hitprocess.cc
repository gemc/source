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
	cndc.date       = "2017-07-13";
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
	sprintf(cndc.database,"/calibration/cnd/Status_LR:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
		cndc.status_L[isec-1][ilay-1][istr-1]=data[row][3];
		cndc.status_R[isec-1][ilay-1][istr-1]=data[row][4];
	}
	
	cout<<"CND:Getting TDC slope"<<endl;
	sprintf(cndc.database,"/calibration/cnd/TDC_conv:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
		cndc.slope_L[isec-1][ilay-1][istr-1]=data[row][3];
		cndc.slope_R[isec-1][ilay-1][istr-1]=data[row][5];
	}
	
	cout<<"CND:Getting attenuation"<<endl;
	sprintf(cndc.database,"/calibration/cnd/Attenuation:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
		cndc.attlen_L[isec-1][ilay-1][istr-1]=data[row][3];
		cndc.attlen_R[isec-1][ilay-1][istr-1]=data[row][5];
	}
	
	cout<<"CND:Getting effective_velocity"<<endl;
	sprintf(cndc.database,"/calibration/cnd/EffV:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
		cndc.veff_L[isec-1][ilay-1][istr-1]=data[row][3];
		cndc.veff_R[isec-1][ilay-1][istr-1]=data[row][5];
	}
	
	cout<<"CND:Getting energy calibration"<<endl;
	sprintf(cndc.database,"/calibration/cnd/Energy:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
		cndc.mip_dir_L[isec-1][ilay-1][istr-1]=data[row][3];
		cndc.mip_indir_L[isec-1][ilay-1][istr-1]=data[row][5];
		cndc.mip_dir_R[isec-1][ilay-1][istr-1]=data[row][7];
		cndc.mip_indir_R[isec-1][ilay-1][istr-1]=data[row][9];
	}
	
	cout<<"CND:Getting u-turn delay"<<endl;
	sprintf(cndc.database,"/calibration/cnd/UturnTloss:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr = data[row][2];
		cndc.uturn_tloss[isec-1][ilay-1][0]=data[row][3];
	}
	
	cout<<"CND:Getting time offset LR"<<endl;
	sprintf(cndc.database,"/calibration/cnd/TimeOffsets_LR:%d",cndc.runNo);
	data.clear(); calib->GetCalib(data,cndc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
		cndc.time_offset_LR[isec-1][ilay-1][0]=data[row][3];
	}
	
	cout<<"CND:Getting time offset layer"<<endl;
	sprintf(cndc.database,"/calibration/cnd/TimeOffsets_layer:%d",cndc.runNo);
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
	if(aHit->isBackgroundHit == 1) return dgtz;

	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);
	
	int sector  = identity[0].id;
	int layer  = identity[1].id;
	int paddle = identity[2].id;
	
	// Get the paddle length: in CND paddles are along z
	double length = aHit->GetDetector().dimensions[0];     // this is actually the half-length! Units: mm
	
	
	// TDC smearing for the time resolution (needs to know which layer it's in) --
	// factor for direct and indirect signals is tuned on cosmic ray measurements and simulations, so should be hard-coded for now:
	
	double uturn_scale[3];     // to account for nominal energy-loss in uturn, for the 3 layers -- used only for the hard-coded resolutions
	uturn_scale[0] = 0.65;
	uturn_scale[1] = 0.6;
	uturn_scale[2] = 0.5;
	double neigh_scale = uturn_scale[layer-1] * exp(-2*length/cm/150.);   // energy scaling due to nominal propagation in neighbour
	
	double sigmaTD = 0.24;                          // direct signal
	double sigmaTN = sigmaTD / sqrt(neigh_scale);   // indirect signal in neighbour
	
	
	// To calculate ADC values:
	
	double dEdxMIP = 1.956;                       // energy deposited by MIP per cm of scintillator material
	double thickness = 3;                         // thickness of each CND paddle
	
	// estimated yield of photoelectrons at the photocathode of PMT:
	// assumes 10,000 photons/MeV, LG length 1.4m with attenuation length 9.5m, 30% losses at junctions and PMT QE = 0.2
	
	double pmtPEYldD = 1210;                      // for the direct signal
	double pmtPEYldN = pmtPEYldD * neigh_scale;   // for the indirect signal, additional loss in u-turn and attenuation in neighbour
	
	
	// Get info about detector material to eveluate Birks effect
	double birks_constant=aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
	
	vector<G4ThreeVector> Lpos = aHit->GetLPos();   // local position wrt centre of the detector piece (ie: paddle): in mm
	vector<double>      Edep = aHit->GetEdep();     // deposited energy in the hit, in MeV
	vector<int> charge = aHit->GetCharges();        // charge for each step
	vector<double> times = aHit->GetTime();
	vector<double> dx = aHit->GetDx();              // step length
	
	unsigned nsteps = times.size();                 // total number of steps in the hit
	
	
	// variables for each step:
	double dUp = 0.;       // distance travelled by the light along the paddle UPSTREAM (direct light)
	double dDown = 0.;     // distance travelled by the light along the paddle DOWNSTREAM (it will then go through the u-turn and along the neighbouring paddle)
	double Edep_B = 0.;    // Energy deposit, scaled in accordance with Birk's effect
	
	double e_up = 0.;      // attenuated energy as it arrives at the upstream edge of the hit paddle
	double e_down = 0.;    // attenuated energy as it arrives at the downstream edge of the hit paddle
	
	
	// Variables for each hit:
	//  int flag_counted = 0;   // Flag to specify that the time of hit has already been determined for this hit (1: yes, 0: no)
	
	double et_D = 0.;         // total energy*time for each hit propagated to upstream end of hit paddle
	double et_N = 0.;         // total energy*time for each hit propagated to upstream end of the neighbour paddle
	
	double etotUp = 0.;       // total energy of hit propagated to the upstream end of paddle
	double etotDown = 0.;     // total energy of hit propagated to the downstream end of the paddle
	double timeD = 0.;        // hit times measured at the upstream edges of the two paddles
	double timeN = 0.;
	
	int TDCmax = 16384;       // max value of the ADC readout
	
	int ADCD = 0;
	int ADCN = 0;
	int TDCD = TDCmax;
	int TDCN = TDCmax;
	
	int ADCL = 0;
	int ADCR = 0;
	int TDCL = TDCmax;
	int TDCR = TDCmax;
	
	
	// Variables from the database:
	
	double attlength_D = 0.;
//	double attlength_N = 0.;
	double v_eff_D = 0.;
	double v_eff_N = 0.;
	double slope_D = 0.;
	double slope_N = 0.;
	double adc_mip_D = 0.;
	double adc_mip_N = 0.;
	int status_D = 0;
	int status_N = 0;
	
	double t_u = cndc.uturn_tloss[sector-1][layer-1][0];
	double t_offset_LR = cndc.time_offset_LR[sector-1][layer-1][0];
	double t_offset_layer = cndc.time_offset_layer[sector-1][layer-1][0];
	
	if (paddle == 1){   // hit is in paddle L
		
		attlength_D = cndc.attlen_L[sector-1][layer-1][0];
//		attlength_N = cndc.attlen_R[sector-1][layer-1][0];
		
		v_eff_D = cndc.veff_L[sector-1][layer-1][0];
		v_eff_N = cndc.veff_R[sector-1][layer-1][0];
		
		slope_D = cndc.slope_L[sector-1][layer-1][0];
		slope_N = cndc.slope_R[sector-1][layer-1][0];
		
		adc_mip_D = cndc.mip_dir_L[sector-1][layer-1][0];
		adc_mip_N = cndc.mip_indir_L[sector-1][layer-1][0];
		
		status_D = cndc.status_L[sector-1][layer-1][0];
		status_N = cndc.status_R[sector-1][layer-1][0];
	}
	
	else if (paddle == 2) {   // hit is in paddle R
		
		attlength_D = cndc.attlen_R[sector-1][layer-1][0];
//		attlength_N = cndc.attlen_L[sector-1][layer-1][0];
		
		v_eff_D = cndc.veff_R[sector-1][layer-1][0];
		v_eff_N = cndc.veff_L[sector-1][layer-1][0];
		
		slope_D = cndc.slope_R[sector-1][layer-1][0];
		slope_N = cndc.slope_L[sector-1][layer-1][0];
		
		adc_mip_D = cndc.mip_dir_R[sector-1][layer-1][0];
		adc_mip_N = cndc.mip_indir_R[sector-1][layer-1][0];
		
		status_D = cndc.status_R[sector-1][layer-1][0];
		status_N = cndc.status_L[sector-1][layer-1][0];
	}
	else cout <<"/n Help, do not recognise paddle number " << paddle << "!!!" << endl;
	
	
	if(tInfos.eTot>0)
	{
		for(unsigned int s=0; s<nsteps; s++)
		{
			// Distances travelled through the paddles to the upstream (dDir) and downstream (dNeigh) edges:
			dUp   = (length + Lpos[s].z());
			dDown = (length - Lpos[s].z());
			//cout<<"Lpos "<<Lpos[s].z()/cm<<endl;
			//cout<<"dUp "<<dUp/cm<<endl;
			//cout<<"dDown "<<dDown/cm<<endl;
			
			// apply Birks effect:
			Edep_B = BirksAttenuation(Edep[s],dx[s],charge[s],birks_constant);
			
			// Calculate attenuated energy which will reach the upstream and downstream edges of the hit paddle:
			e_up   = Edep_B * 0.5 * exp(-dUp/cm/attlength_D);
			e_down = Edep_B * 0.5 * exp(-dDown/cm/attlength_D);
			
			// Integrate energy over entire hit. These values are used for time-smearing:
			etotUp = etotUp + e_up;
			etotDown = etotDown + e_down;
			
			//cout << "step: " << s << " etotUp, etotDown: " << etotUp << ", " << etotDown << endl;
			
			
			/****** Time of hit calculation *******/
			// In all the methods below, the assumption is that the time taken to travel along the light-guides,
			// through PMTs and the electronics to the ADC units is the same for all paddles, therefore this time offset is ignored.
			// Pick whichever method you like:
			
			// Method 1: this takes average time of all the steps with energy deposit above 0:
			//if (Edep[s] > 0.){
			//   timeD = timeD + (times[s] + dUp/cm/v_eff_D) / nsteps;
			//   timeN = timeN + (times[s] + dDown/cm/v_eff_D + t_u + 2*length/cm/v_eff_N) / nsteps;
			//}
			
			// Method 2: This takes the time of the first step (in order of creation) with energy deposit above 0 (can set another threshold):
			// if (flag_counted == 0 && Edep[s] > 0.){
			//   timeD = times[s] + dUp/cm/v_eff_D;
			//   timeN = times[s] + dDown/cm/v_eff_D + t_u + 2*length/cm/v_eff_N;
			//   flag_counted = 1;   // so that subsequent steps are not counted in the hit
			// }
			
			// Method 3: This calculates the total energy * time value at the upstream edges of the hit paddle and its neighbour,
			// will be used to get the energy-weighted average times (should correspond roughly to the peak energy deposit):
			
			// cout<<"times[s] (ns) "<<times[s]/ns<<endl;
			
			et_D = et_D + ((times[s] + dUp/cm/v_eff_D) * e_up);
			et_N = et_N + ((times[s] + dDown/cm/v_eff_D + t_u + 2*length/cm/v_eff_N) * e_down);
		}   // close loop over steps s
		
		/**** The following calculates the time based on energy-weighted average of all step times ****/
		
		timeD = et_D / etotUp;      // sum(energy*time) /  sum(energy)
		timeN = et_N / etotDown;
		
		// "Actual" GEMC hit timings, propagated to paddle edges and not:
		//timeD = tInfos.time + dUp/cm/v_eff_D;
		//timeN = tInfos.time + dDown/cm/v_eff_D + t_u + 2*length/cm/v_eff_N;
		//timeD = tInfos.time;
		//timeN = tInfos.time + t_u + 2*length/cm/v_eff_N;
		
		//	cout<<"timeD "<<timeD<<endl;
		//	cout<<"timeN "<<timeN<<endl;
		//	cout<<"timeD-timeN "<<timeD-timeN<<endl;
		
		// Apply offsets between left and right paddles (these are only applied on paddle 2, which is R):
		
		if (paddle == 2) timeD = timeD + t_offset_LR;
		else if (paddle == 1) timeN = timeN + t_offset_LR;
		
		// Apply global offsets for each paddle-pair (a.k.a. component):
		
		timeD = timeD + t_offset_layer;
		timeN = timeN + t_offset_layer;
		
		/******** end timing determination ***********/
		
		// cout << "Half-length of paddle(cm): "<<length/cm<<endl;
		// cout << "etotUp, etotDown (in MeV): " << etotUp/MeV << ", " << etotDown/MeV << endl;
		// cout << "timeD, timeN (in ns): " << timeD/ns << ", " << timeN/ns << endl;
		// cout << "Reconstructed time (ns): " << ((timeD + timeN - 2*length/cm/v_eff_N - t_u)/2.)/ns << endl;
		// cout << "Reconstructed z (cm, wrt paddle center): " << (length/cm/v_eff_N + t_u + timeD - timeN)*v_eff_D/2./cm << endl;
		
		// Check the full output:
		// cout << "Total steps in this hit: " << nsteps << endl;
		// for (int s=0; s<nsteps; s++)
		//  {
		//    cout << "\n Edep (in MeV): "  << Edep[s] << " time (in ns): " << times[s] << " Lpos-x (in mm): "
		//          << Lpos[s].x() << " Lpos-y: " << Lpos[s].y() << " Lpos-z: " << Lpos[s].z() << " pos-x (in mm): "
		// 		         << pos[s].x() << " pos-y: "  << pos[s].y()  << " pos-z: " << pos[s].z() << endl;
		//  }
		//  cout << "Total for the hit:" << endl;
		//  cout << "etotUp (in MeV): " << etotUp << " etotDown: " << etotDown  << " timeD (in ns): " << timeD << " timeN: " << timeN << endl;
		
		
		/**** Actual digitisation happens here! *****/
		
		if (etotUp > 0.)
		{
			TDCD = (int) ((G4RandGauss::shoot(timeD,sigmaTD/sqrt(etotUp)))/slope_D);
			double npheD = G4Poisson(etotUp*pmtPEYldD);
			double eneD = npheD/pmtPEYldD;
			ADCD = (int) (eneD*adc_mip_D*2./(dEdxMIP*thickness));
		}
		if (etotDown > 0.)
		{
			TDCN = (int) ((G4RandGauss::shoot(timeN,sigmaTN/sqrt(etotDown)))/slope_N);
			double npheN = G4Poisson(etotDown*pmtPEYldN);
			double eneN = npheN/pmtPEYldN;
			ADCN = (int) (eneN*adc_mip_N*2./(dEdxMIP*thickness));
		}
		
		if (TDCD < 0) TDCD = 0;
		else if (TDCD > TDCmax) TDCD = TDCmax;
		if (TDCN < 0) TDCN = 0;
		else if (TDCN > TDCmax) TDCN = TDCmax;
		
		if (ADCD < 0) ADCD = 0;
		if (ADCN < 0) ADCN = 0;
		
	}  // closes tInfos.eTot>0
	
	
	if(verbosity>4)
	{
		cout <<  hd_msg << " layer: " << layer    << ", paddle: " << paddle  << " x=" << tInfos.x/cm << "cm, y=" << tInfos.y/cm << "cm, z=" << tInfos.z/cm << "cm" << endl;
		cout <<  hd_msg << " Etot=" << tInfos.eTot/MeV     << "MeV, average time=" << tInfos.time  << "ns"  << endl;
		cout <<  hd_msg << " TDCD= " << TDCD     << ", TDCN= " << TDCN    << ", ADCD= " << ADCD << ", ADCN= " << ADCN << endl;
	}
	
	// Status flags
	switch (status_D)
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
			cout << " > Unknown CND status: " << status_D << " for sector " << sector << ",  layer " << layer << ", paddle " << paddle << endl;
	}
	switch (status_N)
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
			cout << " > Unknown CND status: " << status_N << " for sector " << sector << ",  layer " << layer << ", paddle " << 3 - paddle << endl;
	}
	
	if (paddle == 1)
	{
		ADCL = (int) ADCD;
		ADCR = (int) ADCN;
		TDCL = (int) TDCD;
		TDCR = (int) TDCN;
	}
	if (paddle == 2)
	{
		ADCL = (int) ADCN;
		ADCR = (int) ADCD;
		TDCL = (int) TDCN;
		TDCR = (int) TDCD;
	}
	
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["layer"]  = layer;
	dgtz["component"] = 1;
	dgtz["ADCL"]   = (int) ADCL;
	dgtz["ADCR"]   = (int) ADCR;
	dgtz["TDCL"]   = (int) TDCL;
	dgtz["TDCR"]   = (int) TDCR;
	
	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	bool rejectHitConditions = false;
	if(rejectHitConditions) {
		writeHit = false;
	}
	
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

