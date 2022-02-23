// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "band_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussT.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

double BARLENGTHS[]  = {163.7,201.9,51.2,51.2,201.9};

static bandHitConstants initializeBANDHitConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	bandHitConstants bhc;

	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return bhc;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}

	bhc.nsector = 6;
	bhc.nlayer = 6;
	bhc.ncomp = 7;

	// database
	bhc.runNo = runno;

	bhc.date       = "2020-07-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		bhc.connection = (string) getenv("CCDB_CONNECTION");
	else
		bhc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";


	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(bhc.connection));


	vector<vector<double> > data;
	int isector, ilayer, icomp;

	//ADD Statustable in the future, F.H 02/08/2021

	//cout<<"BAND:Getting effective velocities"<<endl;
	sprintf(bhc.database,"/calibration/band/effective_velocity:%d:%s%s",bhc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,bhc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isector    = data[row][0];
		ilayer     = data[row][1];
		icomp	   = data[row][2];
		bhc.eff_vel_tdc [isector-1][ilayer-1][icomp-1] = data[row][3];
		bhc.eff_vel_fadc[isector-1][ilayer-1][icomp-1] = data[row][4];
		//printf("%i \t %i \t %i \t %.2f \n", isector, ilayer, icomp, bhc.eff_vel[isector-1][ilayer-1][icomp-1]);

	}

	//cout<<"BAND:Getting attenuation lengths"<<endl;
	sprintf(bhc.database,"/calibration/band/attenuation_lengths:%d:%s%s",bhc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,bhc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isector    = data[row][0];
		ilayer     = data[row][1];
		icomp	   = data[row][2];
		bhc.atten_len[isector-1][ilayer-1][icomp-1] = data[row][3];
		//printf("%i \t %i \t %i \t %.2f \n", isector, ilayer, icomp, bhc.atten_len[isector-1][ilayer-1][icomp-1]);

	}

	//cout<<"BAND:Getting TDC offsets and resolutions"<<endl;
	sprintf(bhc.database,"/calibration/band/paddle_offsets_tdc:%d:%s%s",bhc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,bhc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isector    = data[row][0];
		ilayer     = data[row][1];
		icomp	   = data[row][2];
		bhc.tdc_offset[isector-1][ilayer-1][icomp-1] = data[row][3];
		bhc.tdc_resolution[isector-1][ilayer-1][icomp-1] = data[row][4];
	}

	// These are not in the CCDB
	// Fill with constant values
	//cout<<"BAND:Getting MeV->ADC conversions"<<endl;
	for(isector = 0; isector < bhc.nsector; isector++) {
		for(ilayer = 0; ilayer < bhc.nlayer; ilayer++) {
			for(icomp = 0; icomp < bhc.ncomp; icomp++) {
				bhc.mev_adc[isector][ilayer][icomp] = 1000.; // channels/MeV
			}
		}
	}

	// setting voltage signal parameters
	bhc.vpar[0] = 20;  // delay, ns
	bhc.vpar[1] = 10;  // rise time, ns
	bhc.vpar[2] = 30;  // fall time, ns
	bhc.vpar[3] = 1;   // amplifier

	return bhc;
}

map<string, double> band_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{

	// This is all steps for one identifier

	//cout << "IN INTEGRATE DGT\n";
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	int sector    = identity[0].id;
	int layer     = identity[1].id;
	int component = identity[2].id;
	// int barID = sector*100 + layer*10 + component;

	// You can either loop over all the steps of the hit, or just take the
	// Edep averaged quantities from the trueInfos object:
	trueInfos tInfos(aHit);

	//double Edep = tInfos.eTot;
	//double time = tInfos.time;
	//double x = tInfos.lx / 10.;
	//double y = tInfos.y / 10.;
	//double z = tInfos.z / 10.;
	//int PID = aHit->GetPID();
	//int Z = aHit->GetCharge();

	//// Based on PID, convert MeV Edep to MeVee
	//double E_MeVee = MeVtoMeVee(PID,Z,Edep);

	// Based on x position and time (ToF), convert to PMT L/R time:

	/*
	if(aHit->isBackgroundHit == 1) {

	// background hit has all the energy in the first step. Time is also first step
	double totEdep  = aHit->GetEdep()[0];
	double stepTime = aHit->GetTime()[0];

	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = isector;
	dgtz["layer"]     = ilayer;
	dgtz["component"] = icomponent;
	dgtz["adc"]       = (int) (charge/bhc.fadc_LSB);
	dgtz["tdc"]       = (int) (stepTime*bhc.time_to_tdc);;

	return dgtz;
	}
	*/

	// For each step do birk's laws
	double birks_constant=aHit->GetDetector().GetLogical()->GetMaterial()->GetIonisation()->GetBirksConstant();
	birks_constant = 0.126; // mm/MeV

	vector<G4ThreeVector> 	Lpos 	= aHit->GetLPos();   	// local position wrt centre of the detector piece (ie: paddle): in mm
	vector<G4ThreeVector> 	pos 	= aHit->GetPos();   	// global position, in mm
	vector<double>      	Edep 	= aHit->GetEdep();     	// deposited energy in the hit, in MeV
	vector<int> 		charge	= aHit->GetCharges();   // charge for each step
	vector<int>		pid	= aHit->GetPIDs();	// PIDs for each step
	vector<double> 		times 	= aHit->GetTime();
	vector<double> 		dx 	= aHit->GetDx();        // step length
	unsigned 		nsteps	= times.size();         // total number of steps in the hit

	double 	L 	= BARLENGTHS[sector-1];				 // length of the bar [cm]
	double 	attenL  = bhc.atten_len[sector-1][layer-1][component-1]; // attenuation length of the bar [cm]
	double 	vEff_fadc= bhc.eff_vel_fadc[sector-1][layer-1][component-1]; 	 // effective velocity of bar [cm/ns]
	double 	vEff_tdc = bhc.eff_vel_tdc [sector-1][layer-1][component-1]; 	 // effective velocity of bar [cm/ns]

	//double tL = time + (L/2.-x)/vEff;
	//double tR = time + (L/2.+x)/vEff;

	//cout << time << " " << x << " " << tInfos.x << " " << vEff << "\n";


	double eTotL = 0;
	double eTotR = 0;
	double tL_tdc = 0;
	double tR_tdc = 0;
	double tL_fadc = 0;
	double tR_fadc = 0;
	double xHit = 0;
	double yHit = 0;
	double zHit = 0;
	double rawEtot = 0;
	if( tInfos.eTot > 0 ){
		double et_L_tdc = 0.; // energy-weighted timeL
		double et_R_tdc = 0.; // energy-weighted timeR
		double et_L_fadc = 0.; // energy-weighted timeL
		double et_R_fadc = 0.; // energy-weighted timeR

		double et_X = 0.; // energy-weighted X
		double et_Y = 0.; // energy-weighted Y
		double et_Z = 0.; // energy-weighted Z
		for(unsigned int s=0; s<nsteps; s++){
			// apply Birks effect:
			double Edep_B = BirksAttenuation(Edep[s],dx[s],charge[s],birks_constant);
			//Edep_B = MeVtoMeVee(pid[s],charge[s],Edep[s]);
			rawEtot = rawEtot + Edep[s];

			// Calculate attenuated energy which will reach the upstream and downstream edges of the hit paddle:
			double dL    = (L/2. + Lpos[s].x()/cm);
			double dR    = (L/2. - Lpos[s].x()/cm);
			double e_L   = Edep_B * exp( -dL / attenL);
			double e_R   = Edep_B * exp( -dR / attenL);

			//cout << "step: " << s << "\n";
			//cout << "\t" << Edep[s] << " " << dx[s] << " " << charge[s] << " " << pid[s] << " " << birks_constant << "\n";
			//cout << "\t" << Edep_B << "\n";
			//cout << "\t" << dL << " " << dR << " " << attenL << " " << e_L << " " << e_R << "\n";

			// Integrate energy over entire hit. These values are used for time-smearing:
			eTotL = eTotL + e_L;
			eTotR = eTotR + e_R;

			// Light-output weight the times:
			et_L_tdc = et_L_tdc + (times[s] + dL/vEff_tdc)*e_L;
			et_R_tdc = et_R_tdc + (times[s] + dR/vEff_tdc)*e_R;
			et_L_fadc = et_L_fadc + (times[s] + dL/vEff_fadc)*e_L;
			et_R_fadc = et_R_fadc + (times[s] + dR/vEff_fadc)*e_R;

			et_X = et_X + pos[s].x()/cm * sqrt(e_L*e_R);
			et_Y = et_Y + pos[s].y()/cm * sqrt(e_L*e_R);
			et_Z = et_Z + pos[s].z()/cm * sqrt(e_L*e_R);

		}   // close loop over steps s

		/**** The following calculates the time based on energy-weighted average of all step times ****/

		tL_tdc = et_L_tdc / eTotL;      // sum(energy*time) /  sum(energy)
		tR_tdc = et_R_tdc / eTotR;
		tL_fadc = et_L_fadc / eTotL;      // sum(energy*time) /  sum(energy)
		tR_fadc = et_R_fadc / eTotR;

		xHit = et_X / sqrt(eTotL*eTotR);
		yHit = et_Y / sqrt(eTotL*eTotR);
		zHit = et_Z / sqrt(eTotL*eTotR);
	}



	if( layer == 6 ){ // To avoid infinities due to no effective velocity conversion and no
			  // position determination in veto bars with 1 PMT readout
		tL_tdc = tInfos.time;
		tR_tdc = 0;
		tL_fadc = tInfos.time;
		tR_fadc = 0;
	}

	// Apply simplistic smearing function
	tL_fadc = CLHEP::RandGauss::shoot( tL_fadc, 0.3 );
	tR_fadc = CLHEP::RandGauss::shoot( tR_fadc, 0.3 );
	tL_tdc  = CLHEP::RandGauss::shoot( tL_tdc,  0.3 );
	tR_tdc  = CLHEP::RandGauss::shoot( tR_tdc,  0.3 );


	dgtz["hitn"]      	= (int) hitn;
	dgtz["sector"]    	= (int) sector;
	dgtz["layer"]     	= (int) layer;
	dgtz["component"] 	= (int) component;
	dgtz["ADCL"]		= (int) (1E4 * eTotL);
	//dgtz["amplitudeL"]	= (int) (1E4 * rawEtot);
	dgtz["amplitudeL"]	= (int) (1E4 * xHit);
	dgtz["ADCtimeL"]	= (double) tL_fadc;
	dgtz["TDCL"]		= (int) (1E4 * tL_tdc / 0.02345);
	dgtz["ADCR"]		= (int) (1E4 * eTotR);
	//dgtz["amplitudeR"]	= (int) (1E4 * rawEtot );
	dgtz["amplitudeR"]	= (int) (1E4 * zHit );
	dgtz["ADCtimeR"]	= (double) tR_fadc;
	dgtz["TDCR"]		= (int) (1E4 * tR_tdc / 0.02345);

	//cout << "***************\n";
	//cout << "hitn:\t\t" << hitn << "\n";
	//cout << "sector:\t\t" << sector << "\n";
	//cout << "layer:\t\t" << layer << "\n";
	//cout << "component:\t" << component << "\n";
	//cout << "eTotL:\t\t" << eTotL << "\n";
	//cout << "eTotR:\t\t" << eTotR << "\n";
	//cout << "tL:\t\t" << tL_fadc << "\n";
	//cout << "tR:\t\t" << tR_fadc << "\n";
	//cout << "tL:\t\t" << tL_tdc << "\n";
	//cout << "tR:\t\t" << tR_tdc << "\n";
	//cout << "EffVelTDC:\t" << vEff_tdc << "\n";
	//cout << "EffVelFDC:\t" << vEff_fadc << "\n";
	//cout << "tInfoTime:\t" << tInfos.time << "\n";
	//cout << "x:\t\t" << xHit << "\n";
	//cout << "y:\t\t" << yHit << "\n";
	//cout << "z:\t\t" << zHit << "\n";
	//cout << "***************\n";
	// decide if write an hit or not
	writeHit = true;

	// define conditions to reject hit
	//if(rejectHitConditions) {
	//	writeHit = false;
	//}
	//cout << "RETURNING DGTZ\n";
	return dgtz;
}

vector<identifier>  band_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{

	// This is where it's possible to create two dgtz from 1 hit of geant


	//cout << "IN PROCESS ID\n";
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > band_HitProcess :: chargeTime(MHit* aHit, int hitn)
{

	// All hits within


	//cout << "IN CHARGE TIME\n";
	map< int, vector <double> >  CT;

	vector<double> hitNumbers;
	vector<double> stepIndex;
	vector<double> chargeAtElectronics;
	vector<double> timeAtElectronics;
	vector<double> identifiers;
	vector<double> hardware;
	hitNumbers.push_back(hitn);

	//////////////////////////////////////////////////////////////////
	// COMMENTED BLOCK IS COPIED FROM FT_HODO			//
	// NEED TO IMPLEMENT FOR BAND IF DIGITIZED PULSES ARE REQUIRED	//
	//////////////////////////////////////////////////////////////////

	/*

	// getting identifiers
	vector<identifier> identity = aHit->GetId();


	// use Crystal ID to define IDX and IDY
	int sector    = identity[0].id;
	int layer     = identity[1].id;
	int component = identity[2].id;
	int order = 0; // Always 0

	identifiers.push_back(sector);
	identifiers.push_back(layer);
	identifiers.push_back(component);
	identifiers.push_back(order);

	// getting hardware
	Hardware thisHardware = fthc.TT.getHardware({sector, layer, component, order});
	hardware.push_back(thisHardware.getCrate());
	hardware.push_back(thisHardware.getSlot());
	hardware.push_back(thisHardware.getChannel());

	// Adding pedestal mean and sigma into the hardware as well
	// All of these variables start from 1, therefore -1 is subtracted, e.g. sector-1
	hardware.push_back(fthc.pedestal[sector -1][layer - 1].at(component - 1));
	hardware.push_back(fthc.pedestal_rms[sector -1][layer - 1].at(component - 1));

	trueInfos tInfos(aHit);

	vector<G4ThreeVector> Lpos = aHit->GetLPos();

	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> time = aHit->GetTime();


	for (unsigned int s = 0; s < tInfos.nsteps; s++) {
	// adding shift and spread on time
	double stepTime = time[s] + fthc.time_offset[sector - 1][layer - 1][component - 1] + G4RandGauss::shoot(0., fthc.time_rms[sector - 1][layer - 1][component - 1]);

	// calculate charge and amplitude
	double stepCharge = Edep[s] * fthc.mips_charge[sector - 1][layer - 1][component - 1] / fthc.mips_energy[sector - 1][layer - 1][component - 1];
	double npe_mean = stepCharge / fthc.gain_pc[sector - 1][layer - 1][component - 1];
	double npe = G4Poisson(npe_mean);
	stepCharge = stepCharge * npe / npe_mean;
	//        double amplitude = charge*fthc.gain_mv[isector-1][ilayer-1][icomponent-1]/fthc.gain_pc[isector-1][ilayer-1][icomponent-1];
	//        double fadc      = amplitude/fthc.fadc_LSB;
	double ADC = (stepCharge * fthc.fadc_input_impedence / fthc.fadc_LSB / fthc.ns_per_sample);

	stepIndex.push_back(s);
	chargeAtElectronics.push_back(ADC);
	timeAtElectronics.push_back(stepTime);
	}

*/
	CT[0] = hitNumbers;
	CT[1] = stepIndex;
	CT[2] = chargeAtElectronics;
	CT[3] = timeAtElectronics;
	CT[4] = identifiers;
	CT[5] = hardware;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double band_HitProcess :: voltage(double charge, double time, double forTime)
{
	return PulseShape(forTime, bhc.vpar, charge, time);
}




// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> band_HitProcess :: electronicNoise()
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



map< string, vector <int> >  band_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;

	return MH;
}

void band_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;


	if(bhc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		bhc = initializeBANDHitConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		bhc.runNo = runno;
	}
}


// this static function will be loaded first thing by the executable
bandHitConstants band_HitProcess::bhc = initializeBANDHitConstants(-1);


double band_HitProcess::MeVtoMeVee(int PID, int Z, double E_MeV ){
	double a1, a2, a3, a4;
	if(PID == 11 || PID == -11 || PID == 13 || PID == -13 || PID == 22 || PID == 211 || PID == -211 ){
		// (anti-)electrons, (anti-)muons, gamma, pions
		a1 = 1;
		a2 = 0;
		a3 = 0;
		a4 = 0;
	}
	else if(PID == 2212){ //proton
		a1 = 0.902713 ;
		a2 = 7.55009  ;
		a3 = 0.0990013;
		a4 = 0.736281 ;
	}
	else if(PID == 1000010020){ // deuterium
		a1 = 0.891575 ;
		a2 = 12.2122  ;
		a3 = 0.0702262;
		a4 = 0.782977 ;
	}
	else if(PID == 1000010030){ // tritium
		a1 = 0.881489 ;
		a2 = 15.9064  ;
		a3 = 0.0564987;
		a4 = 0.811916 ;
	}
	else if(PID == 1000020030){ // helium-3
		a1 = 0.803919 ;
		a2 = 34.4153  ;
		a3 = 0.0254322;
		a4 = 0.894859 ;
	}
	else if( Z == 2 ){
		a1 = 0.781501 ;
		a2 = 39.3133  ;
		a3 = 0.0217115;
		a4 = 0.910333 ;
	}
	else if( Z == 3 ){
		a1 = 0.613491 ;
		a2 = 57.1372  ;
		a3 = 0.0115948;
		a4 = 0.951875 ;
	}
	else if( Z == 4 ){
		a1 = 0.435772 ;
		a2 = 45.538   ;
		a3 = 0.0104221;
		a4 = 0.916373 ;
	}
	else if( Z == 5 ){
		a1 = 0.350273 ;
		a2 = 34.4664  ;
		a3 = 0.0112395;
		a4 = 0.912711 ;
	}
	else{ // Treat as C12 and call it good
		a1 = 0.298394 ;
		a2 = 25.5679  ;
		a3 = 0.0130345;
		a4 = 0.908512 ;
	}
	double E_MeVee = a1 * E_MeV  -  a2 * (1-exp(-a3*pow(E_MeV, a4)));
	if( E_MeV < 0. || E_MeVee < 0. ) return 0;
	return E_MeVee;
}

double band_HitProcess::BirksAttenuation(double destep, double stepl, int charge, double birks){

	// Birk's law:
	//  dL = S * dE / (1 + kB * dE/dx + kB_HO * dE/dx * dE/dx )
	//
	//  kB = (1.26-2.02)E-2 [g/(MeVcm2)] for polyvinyltoluene scintillators

	// leading term:
	double kB = 1.26E-2; 	// [g/(MeVcm2])
	double density = 1.023; // [g/cm3]
	kB = kB/density; 	// [cm/MeV]
	// higher order term:
	double kB_HO = 9.6E-6 / (density*density);	// [cm/MeV]

	double response = destep;
	if (birks*destep*stepl*charge != 0.){
		// correction for high charge states
		if( fabs(charge) > 1 ) kB *= (7.2 / 12.6);

		double dedx_cm = (destep) / (stepl/cm); // [MeV/cm]
		response = destep / (1. + kB * dedx_cm + kB_HO * dedx_cm * dedx_cm );
	}
	return response;
}
