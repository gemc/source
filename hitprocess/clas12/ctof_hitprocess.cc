
// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "ctof_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

static ctofConstants initializeCTOFConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false) {
	ctofConstants ctc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return ctc;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	cout << "Entering initializeCTOF" << endl;
	
	ctc.runNo = runno;
	if (getenv("CCDB_CONNECTION") != nullptr)
		ctc.connection = (string) getenv("CCDB_CONNECTION");
	else
		ctc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	ctc.npaddles = 48;
	ctc.thick = 3.0;
	
	ctc.dEdxMIP = 1.956; // muons in polyvinyltoluene
	ctc.dEMIP = ctc.thick * ctc.dEdxMIP;
	ctc.pmtPEYld = 500;
	//	ctc.tdcLSB        = 41.6667; // counts per ns (24 ps LSB)
	
	
	int isec, ilay;
	//	int isec, ilay, istr;
	
	vector<vector<double> > data;
	
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ctc.connection));
	cout << "Connecting to " << ctc.connection << "/calibration/ctof" << endl;
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/attenuation:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting attenuation" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.attlen[isec - 1][ilay - 1][0].push_back(data[row][3]);
		ctc.attlen[isec - 1][ilay - 1][1].push_back(data[row][4]);
	}
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/effective_velocity:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting effective_velocity" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.veff[isec - 1][ilay - 1][0].push_back(data[row][3]);
		ctc.veff[isec - 1][ilay - 1][1].push_back(data[row][4]);
	}
	
	if(accountForHardwareStatus) {
		snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
		cout << "CTOF:Getting status" << endl;
		data.clear();
		calib->GetCalib(data, ctc.database);
		for (unsigned row = 0; row < data.size(); row++) {
			isec = data[row][0];
			ilay = data[row][1];
			//        istr = data[row][2];
			ctc.status[isec - 1][ilay - 1][0].push_back(data[row][3]);
			ctc.status[isec - 1][ilay - 1][1].push_back(data[row][4]);
		}
	}
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/threshold:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting threshold" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.threshold[isec - 1][ilay - 1][0].push_back(data[row][3]);
		ctc.threshold[isec - 1][ilay - 1][1].push_back(data[row][4]);
	}
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/efficiency:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting efficiency" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.efficiency[isec - 1][ilay - 1][0].push_back(data[row][3]);
		ctc.efficiency[isec - 1][ilay - 1][1].push_back(data[row][3]);
	}
	
	
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/gain_balance:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting gain_balance" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.countsForMIP[isec - 1][ilay - 1][0].push_back(data[row][3]);
		ctc.countsForMIP[isec - 1][ilay - 1][1].push_back(data[row][4]);
	}
	
	/* For future use in HitProcess
	 cout<<"Getting time_walk"<<endl;
	 snprintf(ctc.database, sizeof(ctc.database), "/calibration/ctof/time_walk:%d:%s%s",ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	 data.clear() ; calib->GetCalib(data,ctc.database);
	 for(unsigned row = 0; row < data.size(); row++)
	 {
	 isec   = data[row][0]; ilay   = data[row][1]; istr   = data[row][2];
	 ctc.twlk[isec-1][ilay-1][0].push_back(data[row][3]);
	 ctc.twlk[isec-1][ilay-1][1].push_back(data[row][4]);
	 ctc.twlk[isec-1][ilay-1][2].push_back(data[row][5]);
	 ctc.twlk[isec-1][ilay-1][3].push_back(data[row][6]);
	 ctc.twlk[isec-1][ilay-1][4].push_back(data[row][7]);
	 ctc.twlk[isec-1][ilay-1][5].push_back(data[row][8]);
	 }
	 */
	
    snprintf(ctc.database, sizeof(ctc.database), "/calibration/ctof/time_offsets:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting time_offsets" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.toff_UD[isec - 1][ilay - 1].push_back(data[row][3]);
		ctc.toff_RFpad[isec-1][ilay-1].push_back(data[row][4]);
		ctc.toff_P2P[isec-1][ilay-1].push_back(data[row][5]);
	}
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/tdc_conv:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting tdc_conv" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		//        istr = data[row][2];
		ctc.tdcconv[isec - 1][ilay - 1][0].push_back(data[row][3]);
		ctc.tdcconv[isec - 1][ilay - 1][1].push_back(data[row][4]);
	}
	
	
	snprintf(ctc.database, sizeof(ctc.database),  "/geometry/ctof/ctof:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Getting geometry" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		isec = data[row][0];
		ilay = data[row][1];
		double length = data[row][7]/2*cm;
		double offset = data[row][8]*cm;
		ctc.length[isec - 1][ilay - 1].push_back(length);
		ctc.zoffset[isec - 1][ilay - 1].push_back(offset);
		//	cout << data[row][2] << " " << length << " " << offset << endl;
	}
	
	
	snprintf(ctc.database, sizeof(ctc.database),  "/calibration/ctof/tres:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	cout << "CTOF:Setting time resolution" << endl;
	data.clear();
	calib->GetCalib(data, ctc.database);
	for (unsigned row = 0; row < data.size(); row++) {
		double sigma = data[row][3];
		ctc.tres.push_back(sigma);
	}
	
	
	// setting voltage signal parameters
	ctc.vpar[0] = 50; // delay, ns
	ctc.vpar[1] = 10; // rise time, ns
	ctc.vpar[2] = 20; // fall time, ns
	ctc.vpar[3] = 1; // amplifier
	
	// FOR now we will initialize pedestals and sigmas to a random value, in the future
	// they will be initialized from DB
	const double const_ped_value = 101;
	const double const_ped_sigm_value = 2;
	// commands below fill all the elements of ctc.pedestal and ctc.pedestal_sigm with their values (const_ped_value, and const_ped_sigm_value respectively)
	std::fill(&ctc.pedestal[0][0], &ctc.pedestal[0][0] + sizeof (ctc.pedestal) / sizeof (ctc.pedestal[0][0]), const_ped_value);
	std::fill(&ctc.pedestal_sigm[0][0], &ctc.pedestal_sigm[0][0] + sizeof (ctc.pedestal_sigm) / sizeof (ctc.pedestal_sigm[0][0]), const_ped_sigm_value);
	
	// loading translation table
	ctc.TT = TranslationTable("ctofTT");
	
	// loads translation table from CLAS12 Database:
	// Translation table for ctof
	// Crate sector assignments: All FADC readout channels are in the same crate 59, and TDCs are on 60
	// ORDER: 0=Updtream  2=Downstream.
	
	string database = "/daq/tt/ctof:1";
	
	data.clear();
	calib->GetCalib(data, database);
	cout << "  > " << ctc.TT.getName() << " TT Data loaded from CCDB with " << data.size() << " columns." << endl;
	
	// filling translation table
	for (unsigned row = 0; row < data.size(); row++) {
		int crate = data[row][0];
		int slot = data[row][1];
		int channel = data[row][2];
		
		int sector = data[row][3];
		int layer = data[row][4];
		int pmt = data[row][5];
		int order = data[row][6];
		
		// order is important as we could have duplicate entries w/o it
		ctc.TT.addHardwareItem({sector, layer, pmt, order}, Hardware(crate, slot, channel));
	}
	
	
	
	// now connecting to target geometry to get its position
	snprintf(ctc.database, sizeof(ctc.database), "/geometry/target:%d:%s%s", ctc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ctc.database);
	ctc.targetZPos = data[0][3]*cm;
	
	return ctc;
}


map<string, double> ctof_HitProcess::integrateDgt(MHit* aHit, int hitn)
{
	
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	// hit ids
	int sector = 1;
	int layer  = 1;
	int paddle = identity[0].id;
	int side   = identity[1].id;
	
	// TDC conversion factors
	double tdcconv = ctc.tdcconv[sector - 1][layer - 1][side][paddle - 1];
	
	if(aHit->isBackgroundHit == 1) {
		
		// background hit has all the energy in the first step. Time is also first step
		double totEdep = aHit->GetEdep()[0];
		double stepTime = aHit->GetTime()[0];
		
		double adc  = totEdep * ctc.countsForMIP[sector - 1][layer - 1][0][paddle - 1] / ctc.dEMIP ; // no gain as that comes from data already
		double tdc = stepTime/tdcconv;
		
		dgtz["sector"]    = sector;
		dgtz["layer"]     = layer;
		dgtz["component"] = paddle;
		dgtz["ADC_order"] = side;
		dgtz["ADC_ADC"]   = (int) adc;
		dgtz["ADC_time"]  = (tdc*24.0/1000);
		dgtz["ADC_ped"]   = 0;
		
		dgtz["TDC_order"] = side + 2;
		dgtz["TDC_TDC"]   = (int) tdc;
		
		return dgtz;
	}
	
	double length = ctc.length[sector - 1][layer - 1].at(paddle - 1);
	double offset = ctc.zoffset[sector - 1][layer - 1].at(paddle - 1);
	
	
	trueInfos tInfos(aHit);
	
	// Distances from upstream, downstream
	// ctof paddle center is offsetby ctc.zoffset from the CLAS12 target position,
	// so need to  z is also the local coordinate
	// side = 0 or 1,
	double d = length + (1. - 2. * side)*(tInfos.z - offset - ctc.targetZPos); // The distance between the hit and PMT?
	
	// attenuation length
	double attlen = ctc.attlen[sector - 1][layer - 1][side][paddle - 1];
	
	
	double attlen_otherside = ctc.attlen[sector - 1][layer - 1][1 - side].at(paddle - 1);
	
	
	// attenuation factor
	double att = exp(-d / cm / attlen);
	
	// Gain factors to simulate CTOF PMT gain matching algorithm.
	// Each U,D PMT pair has HV adjusted so geometeric mean sqrt(U*D)
	// is independent of counter length, which compensates for
	// the factor exp(-L/2/attlen) where L=full length of bar.
	
	// The old code for the gain is expressed through these two lines,
	//   double gainUp = sqrt(attUp*attDn);    // Line1
	//   double gainDn = gainUp;               // Line2
	// If we look into the expression of the attenuation, the we can quickly
	// find that the gain will be as sqrt(exp(-2*length/cm/attlen))
	
	//double gain = sqrt(exp(-2*length/cm/attlen));
	
	double gain = sqrt(exp(-d / cm / attlen) * exp(-(2 * length - d) / cm / attlen_otherside));
	
	// Attenuated light at PMT
	double energyDepositedAttenuated = tInfos.eTot*att;
	
	double adc = 0.;
	double tdc = 0.;
	// not used anymore
	// double adcu = 0.;
	// double tdcu = 0.;
	// if (ene > 0)
	//	adcu = energyDepositedAttenuated * ctc.countsForMIP[sector - 1][layer - 1][side][paddle - 1] / ctc.dEMIP / gain;
	
	// Fluctuate the light measured by the PMT with
	// Poisson distribution for emitted photoelectrons
	// Treat Up and Dn separately, in case nphe=0
	
	double nphe = G4Poisson(energyDepositedAttenuated * ctc.pmtPEYld);
	double ene = nphe / ctc.pmtPEYld;
	
	if (ene > 0) {
		adc = ene * ctc.countsForMIP[sector - 1][layer - 1][0][paddle - 1] / ctc.dEMIP / gain;
		
		//double            A = ctc.twlk[sector-1][layer-1][0][paddle-1];
		//double            B = ctc.twlk[sector-1][layer-1][1][paddle-1];
		//double            C = ctc.twlk[sector-1][layer-1][2][paddle-1];
		//double   timeWalkUp = A/(B+C*sqrt(adcu));
		double timeWalk = 0.;
		
		double tU = tInfos.time + d/ctc.veff[sector-1][layer-1][side][paddle-1]/cm +  (1. - 2. * side)*ctc.toff_UD[sector-1][layer-1][paddle-1]/2.
		- ctc.toff_RFpad[sector-1][layer-1][paddle-1]
		- ctc.toff_P2P[sector-1][layer-1][paddle-1]
		+ timeWalk;
		
		double t = G4RandGauss::shoot(tU, sqrt(2) * ctc.tres[paddle - 1]);
		// tdcu = tU / tdcconv;
		tdc = t / tdcconv;
	}
	
	if(accountForHardwareStatus) {
		
		// Status flags
		switch (ctc.status[sector - 1][layer - 1][side][paddle - 1]) {
			case 0:
				break;
			case 1:
				adc = 0;
				break;
			case 2:
				tdc = 0;
				break;
			case 3:
				adc = tdc = 0;
				break;
				
			case 5:
				break;
				
			default:
				cout << " > Unknown CTOF status: " << ctc.status[sector - 1][layer - 1][side][paddle - 1] << " for sector " << sector <<
				",  layer " << layer << ", paddle " << paddle << " Side " << side << endl;
		}
	}
	
	int fadc_time = convert_to_precision(tdc*tdcconv);

	
	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = 1;
	dgtz["layer"]     = 1;
	dgtz["component"] = paddle;
	dgtz["ADC_order"] = side;
	dgtz["ADC_ADC"]   = (int) adc;
	dgtz["ADC_time"]  = fadc_time;
	dgtz["ADC_ped"]   = 0;
	
	dgtz["TDC_order"] = side + 2;
	dgtz["TDC_TDC"]   = (int) tdc;
	
	// reject hit if below threshold or efficiency
	if ( energyDepositedAttenuated < ctc.threshold[sector - 1][layer - 1][side][paddle - 1] && applyThresholds ) {
		rejectHitConditions = true;
	}
	
	double random = G4UniformRand();
	if ( random > ctc.efficiency[sector - 1][layer - 1][side][paddle - 1] && applyInefficiencies) {
		rejectHitConditions = true;
	}
	
	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}

vector<identifier> ctof_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector) {	
	
	vector<identifier> yid = id;
	
	yid[0].id_sharing = 1; // Paddle number
	yid[1].id_sharing = 1; // Side: This shows the PMT, whether the upstream one, or the downstream. 0 in the geometry by default
	
	// Check if in the geometry the yid[1].id is 0, this should come from the geometry.
	if (yid[1].id != 0) {
		cout << "*****WARNING***** in ctof_HitProcess :: processID, identifier PM! of the original hit should be 0 " << endl;
		cout << "yid[1].id = " << yid[1].id << endl;
	}
	
	// Now we want to have similar identifiers, but the only difference be id PMT to be 1, instead of 0
	identifier this_id = yid[0];
	yid.push_back(this_id);
	
	// for the second hit, the energy declared is the same but the side id is 1 instead of 0
	this_id = yid[1];
	this_id.id = 1;
	yid.push_back(this_id);
	
	// cout<<"In the ctof processID is.size() = "<<id.size()<<endl;
	// cout<<"In the ctof Id[0] = "<<id[0].id<<endl;
	return yid;
}

// - electronicNoise: returns a vector of hits generated / by electronics.

vector<MHit*> ctof_HitProcess::electronicNoise() {
	vector<MHit*> noiseHits;
	
	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);
	
	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)
	
	//	for(unsigned int p=1; p<10; p++)
	//	{
	//		vector<identifier> thisID;
	//
	//		// for paddle, identifier is only 1 dimensional: paddle ID
	//		identifier thisIdentifier;
	//		thisIdentifier.id = p;
	////		identifier.name = "ctofNoise";
	//
	//		thisID.push_back(thisIdentifier);
	//
	//		MHit* thisNoiseHit = new MHit(10.0*p, (double) p, thisID, p);
	//
	//		noiseHits.push_back(thisNoiseHit);
	//	}
	//
	return noiseHits;
}

map< string, vector <int> > ctof_HitProcess::multiDgt(MHit* aHit, int hitn) {
	map< string, vector <int> > MH;
	
	return MH;
}


// - charge: returns charge/time digitized information / step

map< int, vector <double> > ctof_HitProcess::chargeTime(MHit* aHit, int hitn) {
	map< int, vector <double> > CT;
	
	vector<double> hitNumbers;
	vector<double> stepIndex;
	vector<double> chargeAtElectronics;
	vector<double> timeAtElectronics;
	vector<double> identifiers;
	vector<double> hardware;
	hitNumbers.push_back(hitn);
	
	vector<identifier> identity = aHit->GetId();
	
	int sector = 1;
	int layer = 1;
	int paddle = identity[0].id;
	int side = identity[1].id;
	
	identifiers.push_back(sector);
	identifiers.push_back(layer);
	identifiers.push_back(paddle);
	identifiers.push_back(side);
	
	// getting hardware
	Hardware thisHardware = ctc.TT.getHardware({sector, layer, paddle, side});
	hardware.push_back(thisHardware.getCrate());
	hardware.push_back(thisHardware.getSlot());
	hardware.push_back(thisHardware.getChannel());
	
	
	// Adding pedestal mean and sigma into the hardware as well
	// All of these variables start from 1, therefore -1 is subtracted, e.g. sector-1
	hardware.push_back(ctc.pedestal[paddle - 1][side]);
	hardware.push_back(ctc.pedestal_sigm[paddle - 1][side]);
	
	//==== Testing Rafo ====
	//        cout<<"identity[0].id = "<<identity[0].id<<endl;
	//        cout<<"identity[1].id = "<<identity[1].id<<endl;
	
	double length = ctc.length[sector - 1][layer - 1].at(paddle - 1);
	double offset = ctc.zoffset[sector - 1][layer - 1].at(paddle - 1);
	
	trueInfos tInfos(aHit);
	
	// attenuation length
	double attlen = ctc.attlen[sector - 1][layer - 1][side].at(paddle - 1);
	
	double attlen_otherside = ctc.attlen[sector - 1][layer - 1][1 - side].at(paddle - 1);
	
	
	vector<G4ThreeVector> Pos = aHit->GetPos();
	
	// Vector of Edep and time of the hit in each step
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> time = aHit->GetTime();
	
	
	for (unsigned int s = 0; s < tInfos.nsteps; s++) {
		//pmt = 0 or 1,
		
		double d = length + (1. - 2. * side)*(Pos[s].z() - offset); // The distance between the hit and PMT?
		// attenuation factor
		double att = exp(-d / cm / attlen);
		
		// Gain factors to simulate CTOF PMT gain matching algorithm.
		// Each U,D PMT pair has HV adjusted so geometeric mean sqrt(U*D)
		// is independent of counter length, which compensates for
		// the factor exp(-L/2/attlen) where L=full length of bar.
		
		// The old code for the gain is expressed through these two lines,
		//   double gainUp = sqrt(attUp*attDn);    // Line1
		//   double gainDn = gainUp;               // Line2
		// If we look into the expression of the attenuation, the we can quickly
		// find that the gain will be as sqrt(exp(-2*length/cm/attlen))
		
		//double gain = sqrt(exp(-2*length/cm/attlen));
		
		double gain = sqrt(exp(-d / cm / attlen) * exp(-(2 * length - d) / cm / attlen_otherside));
		
		// Attenuated light at PMT
		double ene = tInfos.eTot*att;
		
		// In Mode1 we don't need unsmeared
		
		double nphe = G4Poisson(ene * ctc.pmtPEYld);
		ene = nphe / ctc.pmtPEYld;
		
		if (ene > 0) {
			double adc = ene * ctc.countsForMIP[sector - 1][layer - 1][side][paddle - 1] / ctc.dEMIP / gain;
			//double            A = ctc.twlk[sector-1][layer-1][0][paddle-1];
			//double            B = ctc.twlk[sector-1][layer-1][1][paddle-1];
			//double            C = ctc.twlk[sector-1][layer-1][2][paddle-1];
			//double   timeWalkUp = A/(B+C*sqrt(adcu));
			double timeWalk = 0.;
			
			double tU = time[s] + d/ctc.veff[sector-1][layer-1][side][paddle-1]/cm +  (1. - 2. * side)*ctc.toff_UD[sector-1][layer-1][paddle-1]/2.
			- ctc.toff_RFpad[sector-1][layer-1][paddle-1]
			- ctc.toff_P2P[sector-1][layer-1][paddle-1]
			+ timeWalk;
			
			
			double t = G4RandGauss::shoot(tU, sqrt(2) * ctc.tres[paddle - 1]);
			
			stepIndex.push_back(s); // Since it is going to be only one hit, i.e. only one step
			chargeAtElectronics.push_back(adc);
			timeAtElectronics.push_back(t);
		}
		
	}
	
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

double ctof_HitProcess::voltage(double charge, double time, double forTime) {
	//	return 0.0;
	//return DGauss(forTime, ctc.vpar, charge, time);
	return PulseShape(forTime, ctc.vpar, charge, time);
}

void ctof_HitProcess::initWithRunNumber(int runno) {
	
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if (ctc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ctc = initializeCTOFConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		ctc.runNo = runno;
	}
}


// this static function will be loaded first thing by the executable
ctofConstants ctof_HitProcess::ctc = initializeCTOFConstants(-1);




