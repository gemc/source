#include <math.h>
#include <random>

// G4 Headers
#include "Randomize.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "rtpc_hitprocess.h"

// 2023.10.03 copied from rich_hitprocess.cc. Added in case for the future needed. (y.-c.)
static rtpcConstants initializeRTPCConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	rtpcConstants rtpcc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return rtpcc;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	cout << "Entering initializeRTPCConstants" << endl;
	
	// database
	rtpcc.runNo = runno;

	if(getenv ("CCDB_CONNECTION") != nullptr)
		rtpcc.connection = (string) getenv("CCDB_CONNECTION");
	else
		rtpcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(rtpcc.connection));
	cout << "Connecting to " << rtpcc.connection << "/calibration/rtpc" << endl;
	
	
	// get time parameter from CCDB (instead of hard code in the hitprocess file.)
	vector<vector<double> > data;
	
	cout << "RTPC:Getting drift parameters" << endl;
	snprintf(rtpcc.database, sizeof(rtpcc.database), "/calibration/rtpc/time_parms:%d:%s%s", rtpcc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear();
	calib -> GetCalib(data, rtpcc.database);
	int iComponent = 0;
	for (unsigned row = 0; row < data.size(); row++) {
		//iSector = data[row][0];
		//iLayer = data [row][1];
		iComponent = data[row][2];
		rtpcc.z0[iComponent - 1] = data[row][3];
		rtpcc.z2[iComponent - 1] = data[row][5];
		rtpcc.z4[iComponent - 1] = data[row][7];

		//cout << row << ": " << data[row][0] << "  " << data[row][1] << "  " << data[row][2] << "  " << data[row][3] << "  " << data[row][4] << "  " << data[row][5] << "  " << data[row][6] << "  " << data[row][7] << endl;
	} // row++
	// the last row is used for tdiff (not used in old hitprocess y.-c. Oct30.2023)

/*
	cout << "RTPC:Getting gain balance" << endl; // keep for future or reference
	snprintf(rtpcc.database, sizeof(rtpcc.database), "/calibration/rtpc/gain_balance:%d:%s%s", rtpcc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear();
	calib -> GetCalib(data, rtpcc.database);
	iComponent = 0;
	int iLayer = 0;
	for (unsigned row = 0; row < data.size(); row++) { // data.size() = 180*96 = 17280
		//iSector = data[row][0];
		iLayer = data[row][1];
		iComponent = data[row][2];
		rtpcc.gain_balance[iComponent][iLayer] = data[row][3]; // loop layer first, and then component in gain_balance table
	} // row++
*/

	//rtpcc.date       = "2016-03-15";  // useful?
	//rtpcc.variation  = "default";
	
	// These are not in the CCDB (needed?  y.-c. Oct30.2023)
	
	return rtpcc;
}



map<string, double> rtpc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{

	rejectHitConditions = false;
	writeHit = true;
	
	// Drift time from first GEM to readout pad
	float t_2GEM2 = 0.0; // 169.183;
	float t_2GEM3 = 0.0; //222.415;
	float t_2PAD = 0.0; //414.459;
	
	// parameters to change for gas mixture and potential
	// gas mixture = He:CO2 at 85.5:14.5
	// potential = 4450V
	
	// --------- Updated 25 June 2020 B-field parameters ----------- //
	//double a_t0 = 1470.0; // -0.607056;
	//double a_t1 = 0.009; // -0.00138137;
	//double a_t2 = -9.0e05; // 1.6494e-05; //Gave crazy big -ve values to tdc or time 
	//double a_t2 = -9.0e-05; // 1.6494e-05;
	
	//double b_t0 = 3290.0; // 3183.62;
	//double b_t1 = -0.0138; // -0.0269385;
	//double b_t2 = -0.000332; // -0.00037645;
	
	//double c_t0 = 0.0; // -0.0200932;
	//double c_t1 = 0.0; // 0.00319441;
	//double c_t2 = 0.0; // -8.75299e-06;
	
	//double a_phi0 = -0.113; // 0;
	//double a_phi1 = -1.05e-6; // 0;
	//double a_phi2 = 1.58e-8; // 0;
	
	//double b_phi0 = -0.9697; // 1.22249;
	//double b_phi1 = 9.48e-6; // -2.3708e-05;
	//double b_phi2 = 9.9e-8; // -1.11055e-07;
	
	//double c_phi0 = 0.0; // -1.13197;
	//double c_phi1 = 0.0; // 0.000118575;
	//double c_phi2 = 0.0; // 2.50094e-08;

	//cout << "In integrateDgt: a_t para = " << a_t0 << "  " << a_t1 << "  " << a_t2 << endl;
	//for (int i = 0; i < 7; i++) {
	//	cout << "In integrateDgt: CCDB cons " << i <<" = " << rtpcc.z0[i] << "  " << rtpcc.z2[i] << "  " << rtpcc.z4[i] << endl;
	//}
	//cout << " ------------------------------------- " << endl;
	
	// Read drift parameters from CCDB -->  checked! Same as the hard code valuse above.
	double a_t0 = rtpcc.z0[0];
	double a_t1 = rtpcc.z2[0];
	double a_t2 = rtpcc.z4[0];

	double b_t0 = rtpcc.z0[1];
	double b_t1 = rtpcc.z2[1];
	double b_t2 = rtpcc.z4[1];
	  // phi offset (delta_phi)
	double a_phi0 = rtpcc.z0[3];
	double a_phi1 = rtpcc.z2[3];
	double a_phi2 = rtpcc.z4[3];
	  // tan_phi
	double b_phi0 = rtpcc.z0[4];
	double b_phi1 = rtpcc.z2[4];
	double b_phi2 = rtpcc.z4[4];

	
	// Diffusion parameters
	diff_at = 388.7449859;
	diff_bt = -4.33E+01;
	
	diff_aphi = 6.00E-06;
	diff_bphi = 2.00E-06;
	
	a_z = 0.035972097;
	b_z = -0.000739386;
	
	float t_2END = t_2GEM2 + t_2GEM3 + t_2PAD;
	sigma_t_gap = sqrt(pow(sigma_t_2GEM2,2)+pow(sigma_t_2GEM3,2)+pow(sigma_t_2PAD,2));
	phi_2END = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
	sigma_phi_gap = sqrt(pow(sigma_phi_2GEM2,2)+pow(sigma_phi_2GEM3,2)+pow(sigma_phi_2PAD,2));
	
	
	TPC_TZERO = 0.0;
	
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
	
	if(tInfos.eTot > 0) {

		// Idea: After collecting the geant4 steps as a hit in processID, re-evaluate the destination on the readout padboard of this hit. Therefore, the final col and row write out in RTPC::adc bank might not be the same as the original(in processID). Futhermore, there might be some different hits have same col and row in RTPC::adc bank.
		// y.-c. comment in Dec. 23. 2023
	
		// total energy deposit for a hit
		double DiffEdep = tInfos.eTot;
		int adc = ((int)100000.0*DiffEdep); //cludge for tiny ADC numbers

		//convert (x0,y0,z0) into (r0,phi0,z0)
		z_cm = tInfos.lz/10.0;  // mm -> cm

		double r0 = (sqrt(tInfos.lx*tInfos.lx + tInfos.ly*tInfos.ly))/10.0;  //in cm
		double phi0_rad = atan2(tInfos.ly, tInfos.lx); //return (-Pi, + Pi)
		if (phi0_rad < 0.) phi0_rad += 2.0*PI; // convert to (0, 2PI)
		//if (phi0_rad >= 2.0*PI) phi0_rad -= 2.0*PI; useless because -pi <= phi0_rad <= pi

		// calculate drift time of a hit
		  // all in cm
		a_t = a_t0 + a_t1*(pow(z_cm,2)) + a_t2*(pow(z_cm,4));  // t_offset
		b_t = b_t0 + b_t1*(pow(z_cm,2)) + b_t2*(pow(z_cm,4));  // T_max
		     
		// --------------------- Addition of Diffusion ----------------- //
		// calculate drift time [ns] to first GEM
		double t_drift = a_t + b_t*(7.0*7.0 - r0*r0)/40.0;
		// determine sigma of drift time [ns]
		double t_diff = sqrt(diff_at*(7.0 - r0) + diff_bt*(7.0 - r0)*(7.0 - r0));

		// find t_s2pad by gaussians
		double t_s2pad = G4RandGauss::shoot(t_drift, t_diff);
		double t_gap = G4RandGauss::shoot(t_2END, sigma_t_gap);
		
		// time shift
		shift_t = timeShift_map.find(otids[0])->second; //To ensure secondaries have the same shift_t as primaries or acencestors
		//tdc = t_s2pad + t_gap + shift_t; //kp:
		double tdc = t_s2pad + t_gap;//+shift_t; scattered electron and spectator proton do not have shift_t



		// calculate drift angle to first GEM at 7 cm [rad]
		  // all in cm 
		a_phi = a_phi0 + a_phi1*(pow(z_cm,2)) + a_phi2*(pow(z_cm,4));  //delta_phi
		b_phi = b_phi0 + b_phi1*(pow(z_cm,2)) + b_phi2*(pow(z_cm,4));  //tan_theta
		  // phi_drift = delta_phi + tanTheta*ln(r_max/r);
		double phi_drift = a_phi + b_phi*log(7.0/r0);

		  // determine sigma of drift angle [rad]
		double phi_diff = sqrt(diff_aphi*(7.0 - r0) + diff_bphi*(7.0 - r0)*(7.0 - r0));
	
		  // find delta_phi by gaussians
		double delta_phi = G4RandGauss::shoot(phi_drift, phi_diff);
		double phi_gap = G4RandGauss::shoot(phi_2END, sigma_phi_gap);

		double phi_rad= phi0_rad + delta_phi + phi_gap;   //phi at pad pcb board
		if (phi_rad < 0.0) phi_rad += 2.0*PI; // return to (0, 2PI)
		if (phi_rad >= 2.0*PI) phi_rad -= 2.0*PI; // return to (0, 2PI)


		// calculate drift in z [cm]
		double z_drift = 0.0;
  
		  // determine sigma in z [cm]
		double z_diff = sqrt(a_z*(7.0 - r0)+b_z*(7.0 - r0)*(7.0 - r0));

		  // find delta_z by gaussians
		double delta_z = G4RandGauss::shoot(z_drift, z_diff);
  
		//double z_pos = (z_cm + delta_z)*10.0; // mm
		double z_pos = z_cm*10.0 + delta_z; // mm



		// convert physical position into row & col of each pad.
		int col = -999;
		int row = -999;

		row = ceil(phi_rad/phi_per_pad);
		float z_shift = (row - 1)%4;
  
		float col_min = -RTPC_L/2.0 + z_shift;
		float col_max = RTPC_L/2.0 + z_shift;
  
		if (z_pos < col_min || z_pos > col_max) {row = -999; col = -999;}
		else {
		  row = ceil(phi_rad/phi_per_pad);
		  col = ceil((z_pos + RTPC_L/2.0 - z_shift)/PAD_L);
		} // consistent with PadVector.java in coatjava
	

//		row = identity[0].id;
//		col = identity[1].id;
		
		//Mohammad: kill all the hits in the last three rows that correspond to the seam effect
		if(row > 177){row = -999; col = -999;}


		// write out to adc bank 
		dgtz["sector"]        = 1;
		dgtz["layer"]         = col;
		dgtz["component"]     = row;
		dgtz["ADC_order"]     = 0;
		dgtz["ADC_time"]      = tdc;
		dgtz["ADC_ADC"]       = (int) adc;
		dgtz["ADC_ped"]       = 0;
		
		//dgtz["ADC_hitn"]      = (int) hitn; // not define in CLARA see in data.json

		//reject hit write out to the adc bank, 30 was a number assigned for test. 
		// apply with the efficiency table, for instance, fot protential future impovement 
		//if (adc < 30)
		//	rejectHitConditions = true;

	} // tInfos.eTot > 0


	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}



vector<identifier>  rtpc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	//cout << " In processID ***************" << endl;
	
	// parameters to change for gas mixture and potential
	// gas mixture = He:CO2 at 85.5:14.5
	// potential = 4450V
	
	// --------- Updated 14 Dec 2020 B-field parameters
	
	// Read drift angle parameters from CCDB -->  checked! Same as the hard code valuse above.
	  // phi offset (delta_phi)
	double a_phi0 = rtpcc.z0[3];
	double a_phi1 = rtpcc.z2[3];
	double a_phi2 = rtpcc.z4[3];
	  // tan_phi
	double b_phi0 = rtpcc.z0[4];
	double b_phi1 = rtpcc.z2[4];
	double b_phi2 = rtpcc.z4[4];

	//cout << "row 4: " << a_phi0 << "  " << a_phi1 << "  " << a_phi2 << endl;
	//cout << "row 5: " << b_phi0 << "  " << b_phi1 << "  " << b_phi2 << endl;

	// Diffusion parameters
	diff_aphi = 6.00E-06;
	diff_bphi = 2.00E-06;
	
	
	a_z = 0.035972097;
	b_z = -0.000739386;

	phi_2END = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
	sigma_phi_gap = sqrt(pow(sigma_phi_2GEM2,2)+pow(sigma_phi_2GEM3,2)+pow(sigma_phi_2PAD,2));
	
	vector<identifier> yid = id;
	
	// Get Coordinates of interaction from Geant4 (aStep)
	G4ThreeVector xyz = aStep -> GetPostStepPoint() -> GetPosition(); // global
	G4ThreeVector Lxyz = aStep -> GetPreStepPoint() -> GetTouchableHandle() -> GetHistory() -> GetTopTransform().TransformPoint(xyz); // local
	
	//int chan = 0;
	double LposX = Lxyz.x();
	double LposY = Lxyz.y();
	double LposZ = Lxyz.z();
	
	z_cm = LposZ/10.0; // convert mm to cm
	
	// all in cm, because the drift parameters were extracted in cm.
	a_phi = a_phi0 + a_phi1*(pow(z_cm,2)) + a_phi2*(pow(z_cm,4)); // delta_phi
	b_phi = b_phi0 + b_phi1*(pow(z_cm,2)) + b_phi2*(pow(z_cm,4)); // tan_theta
	//c_phi=c_phi2*(pow(z_cm,4)) + c_phi1*(pow(z_cm,2)) + c_phi0;
	

	//convert (x0,y0,z0) into (r0,phi0,z0)
	double r0 = (sqrt(LposX*LposX + LposY*LposY))/10.0;  //in cm
	
	double phi0_rad = atan2(LposY,LposX); //return (-Pi, + Pi)
	if (phi0_rad < 0.) phi0_rad += 2.0*PI; // convert to (0, 2PI)
	//if (phi0_rad >= 2.0*PI) phi0_rad -= 2.0*PI; useless because -pi <= phi0_rad <= pi

	
	// -------------------------------- Addition of Diffusion -----------------------------
	
	// calculate drift angle to first GEM at 7 cm [rad]
	  // phi_drift = delta_phi + tanTheta*ln(r_max/r);
	double phi_drift = a_phi + b_phi*log(7.0/r0);
	//double phi_drift = a_phi + b_phi*log(7.0/r0) + c_phi*((1.0/(r0*r0)) - (1.0/(49.0)));
	
	  // determine sigma of drift angle [rad]
	double phi_diff = sqrt(diff_aphi*(7.0 - r0) + diff_bphi*(7.0 - r0)*(7.0 - r0));
	
	// calculate drift in z [cm]
	double z_drift = 0.0;
	
	// determine sigma in z [cm]
	double z_diff = sqrt(a_z*(7.0 - r0)+b_z*(7.0 - r0)*(7.0 - r0));
	
	// find delta_phi by gaussians
	double delta_phi = G4RandGauss::shoot(phi_drift, phi_diff);
	double delta_z = G4RandGauss::shoot(z_drift, z_diff);
	double phi_gap = G4RandGauss::shoot(phi_2END, sigma_phi_gap);
	
	
	// ------------------------------------------------------------------------------------
	
	double phi_rad= phi0_rad + delta_phi + phi_gap;   //phi at pad pcb board
	if (phi_rad < 0.0) phi_rad += 2.0*PI; // return to (0, 2PI)
	if (phi_rad >= 2.0*PI) phi_rad -= 2.0*PI; // return to (0, 2PI)
	
	//double z_pos = (z_cm + delta_z)*10.0; // mm
	double z_pos = z_cm*10.0 + delta_z; // mm
	int col = -999;
	int row = -999;
	
	row = ceil(phi_rad/phi_per_pad);
	float z_shift = (row - 1)%4;
	
	float col_min = -RTPC_L/2.0 + z_shift;
	float col_max = RTPC_L/2.0 + z_shift;
	
	if (z_pos < col_min || z_pos > col_max) {row = -999; col = -999;}
	else {
	  row = ceil(phi_rad/phi_per_pad);
	  col = ceil((z_pos + RTPC_L/2.0 - z_shift)/PAD_L);
	}
	
	yid[0].id = row;
	yid[1].id = col;

/*	cout << " In processID ***************" << endl;
	cout << " input id size: " << id.size() << " name: " << id[0].name << ", " << id[1].name << endl;
	cout << "row: " << row << " | " << col << endl;
	cout << "id: " << yid[0].id << " | " << yid[1].id << " || rule: " << yid[0].rule << " id_sharing: " << yid[0].id_sharing << " | " << yid[1].id_sharing << " time: " <<  yid[0].time << " | " << yid[1].time << " name: " << yid[0].name << ", " << yid[1].name << endl;
*/
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



void rtpc_HitProcess::initWithRunNumber(int runno)
{
	//cout << "*** In initWithRunNumber: runno = " << "   " << runno << endl;
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;

	if(rtpcc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		rtpcc = initializeRTPCConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		rtpcc.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
// setup 11 instead of -1 (in the clas12Tags) because z0, z2, z4 are also used in processID, they need to be loaded first. Moreover, only run 11 in CCDB has appropricate parameters for simulation. 
rtpcConstants rtpc_HitProcess::rtpcc = initializeRTPCConstants(-1);

// add class and member function which can read constants from ccdb in rtpc_hitprocess.X files for the future needed, but comment out them for now (2023.10.03)
// 
