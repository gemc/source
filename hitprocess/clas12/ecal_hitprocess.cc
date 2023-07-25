// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "ecal_hitprocess.h"

static ecConstants initializeECConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	ecConstants ecc;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return ecc;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	int isec,ilay;
	
	// database
	ecc.runNo      = runno;
	
	if(getenv ("CCDB_CONNECTION") != nullptr) {
		ecc.connection = (string) getenv("CCDB_CONNECTION");
	} else {
		ecc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	}
	
	ecc.ADC_GeV_to_evio     = 1./10000.; // MIP based calibration is nominally 10 channels/MeV
	ecc.pmtQE               = 0.27    ;
	ecc.pmtDynodeGain       = 4.0     ;
	
	//  Fluctuations in PMT gain distributed using Gaussian with
	//  sigma = sqrt(npe)/SNR where 1/SNR = sqrt[(1 + 1/(ecc.pmtDynodeGain-1)) npe=number of photoelectrons
	//  Adapted from G-112 (pg. 174) of RCA PMT Handbook.
	
	ecc.pmtFactor           = sqrt(1 + 1/(ecc.pmtDynodeGain-1));	
	
	// The callibration data will be filled in this vector data
	vector<vector<double> > data;
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ecc.connection));
	
	// ======== Initialization of EC gains ===========
	snprintf(ecc.database, sizeof(ecc.database), "/calibration/ec/gain:%d:%s%s", ecc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ecc.database);
	
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec = data[row][0]; ilay = data[row][1];
		ecc.gain[isec-1][ilay-1].push_back(data[row][3]);
	}
	
	// ========= Initializations of attenuation lengths ========
	snprintf(ecc.database, sizeof(ecc.database), "/calibration/ec/attenuation:%d:%s%s", ecc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ecc.database);
	
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec = data[row][0]; ilay = data[row][1];
		ecc.attlen[isec-1][ilay-1][0].push_back(data[row][3]);
		ecc.attlen[isec-1][ilay-1][1].push_back(data[row][5]);
		ecc.attlen[isec-1][ilay-1][2].push_back(data[row][7]);
	}
	
	// ========== Initialization of timings ===========
	snprintf(ecc.database, sizeof(ecc.database), "/calibration/ec/timing:%d:%s%s", ecc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ecc.database);
	
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec = data[row][0]; ilay = data[row][1];
		ecc.timing[isec-1][ilay-1][0].push_back(data[row][3]);
		ecc.timing[isec-1][ilay-1][1].push_back(data[row][4]);
		ecc.timing[isec-1][ilay-1][2].push_back(data[row][5]);
		ecc.timing[isec-1][ilay-1][3].push_back(data[row][6]);
		ecc.timing[isec-1][ilay-1][4].push_back(data[row][7]);
	}
	
	// ========== Initialization of timing offset ===========
	snprintf(ecc.database, sizeof(ecc.database), "/calibration/ec/tdc_global_offset:%d:%s%s", ecc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data, ecc.database);
	ecc.tdc_global_offset = data[0][3];
	
	// ======== Initialization of EC effective velocities ===========
	snprintf(ecc.database, sizeof(ecc.database), "/calibration/ec/effective_velocity:%d:%s%s", ecc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear(); calib->GetCalib(data,ecc.database);
	
	for(unsigned row = 0; row < data.size(); row++)
	{
		isec = data[row][0]; ilay = data[row][1];
		ecc.veff[isec-1][ilay-1].push_back(data[row][3]);
	}
	
	// ======== Initialization of EC status  ===========
	if(accountForHardwareStatus) {
		snprintf(ecc.database, sizeof(ecc.database),  "/calibration/ec/status:%d:%s%s", ecc.runNo, digiVariation.c_str(), timestamp.c_str());
		data.clear();
		calib->GetCalib(data, ecc.database);
		for (unsigned row = 0; row < data.size(); row++) {
			isec = data[row][0]; ilay = data[row][1];
			ecc.status[isec-1][ilay-1].push_back(data[row][3]);
		}	
	}
	
	// =========== Initialization of FADC250 related informations, pedestals, nsa, nsb ======================
	// FOR now we will initialize pedestals and sigmas to a random value,
	// in the future they should come from CCDB
	const double const_ped_value = 101;
	const double const_ped_sigm_value = 2;
	// commands below fill all the elements of ecc.pedestal and ecc.pedestal_sigm with their values (const_ped_value, and const_ped_sigm_value respectively)
	std::fill(&ecc.pedestal[0][0][0], &ecc.pedestal[0][0][0] + sizeof(ecc.pedestal)/sizeof(ecc.pedestal[0][0][0]), const_ped_value);
	std::fill(&ecc.pedestal_sigm[0][0][0], &ecc.pedestal_sigm[0][0][0] + sizeof(ecc.pedestal_sigm)/sizeof(ecc.pedestal_sigm[0][0][0]), const_ped_sigm_value);
	
	// setting voltage signal parameters
	ecc.vpar[0] = 0.;  // delay, ns
	ecc.vpar[1] = 2.8; // rise time, ns
	ecc.vpar[2] = 20;  // fall time, ns
	ecc.vpar[3] = 1;   // amplifier
	
	
	// loading translation table
	ecc.TT = TranslationTable("ecTT");
	
	// loads translation table from CLAS12 Database:
	// Translation table for EC (ECAL+PCAL).
	// Crate sector assignments: ECAL/FADC=1,7,13,19,25,31 ECAL/TDC=2,8,14,20,26,32
	// PCAL/FADC=3,9,15,21,27,33 PCAL/TDC=4,10,16,22,28,34.
	// ORDER: 0=FADC 2=TDC.
	
	string database   = "/daq/tt/ec:1";
	
	data.clear(); calib->GetCalib(data, database);
	cout << "  > " << ecc.TT.getName() << " TT Data loaded from CCDB with " << data.size() << " columns." << endl;
	
	// filling translation table
	for(unsigned row = 0; row < data.size(); row++)
	{
		int crate   = data[row][0];
		int slot    = data[row][1];
		int channel = data[row][2];
		
		int sector  = data[row][3];
		int layer   = data[row][4];
		int pmt     = data[row][5];
		int order   = data[row][6];
		
		// order is important as we could have duplicate entries w/o it
		ecc.TT.addHardwareItem({sector, layer, pmt, order}, Hardware(crate, slot, channel));
	}
	cout << "  > Data loaded in translation table " << ecc.TT.getName() << endl;
	
	return ecc;
}


// Process the ID and hit for the EC using EC scintillator slab geometry instead of individual strips.
map<string, double> ecal_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	int sector = identity[0].id;
	int layer  = identity[1].id; // layer=1-3 (PCAL) 4-9 (ECAL). Layer = view for pcal, ecinner, ecouter
	int strip  = identity[2].id;
	// pcal
	int view   = layer;
	
	bool isPCAL = layer < 4 ;
	
	double time_in_ns = 0;

	
	// layer = 1, 2 stays the same
	// subtract 3 from ec inner
	// subtract 6 from ec outer
	
	if (layer > 3 && layer <7) {
		// ec inner (stack 1)
		view = layer - 3;
	} else if (layer >= 7) {
		// ec outer (stack 2)
		view = layer - 6;
	}
	
	// Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
	// Different for EC and PCAL
	double pmtPEYld = 3.5 ;
	if (isPCAL) {
		// pcal
		pmtPEYld  = 11.5 ;
	}

	double a1   = ecc.timing[sector-1][layer-1][1][strip-1]; // tdc conversion

	if(aHit->isBackgroundHit == 1) {
		
		// background hit has all the energy in the first step. Time is also first step
		double totEdep = aHit->GetEdep()[0];
		double stepTime = aHit->GetTime()[0];
		double adc  = totEdep / ecc.ADC_GeV_to_evio ; // no gain as that comes from data already
		int tdc = stepTime / a1 ;
		
		dgtz["hitn"]      = hitn;
		dgtz["sector"]    = sector;
		dgtz["layer"]     = layer;
		dgtz["component"] = strip;
		dgtz["ADC_order"] = 0;
		dgtz["ADC_ADC"]   = (int) adc;
		dgtz["ADC_time"]  = convert_to_precision(stepTime);
		dgtz["ADC_ped"]   = 0;
		
		dgtz["TDC_order"] = layer < 4 ? 2 : 1;  // 1 ECAL, 2 PCAL
		dgtz["TDC_TDC"]   = tdc;
		
		return dgtz;
	}
	
	HCname = "ECAL Hit Process";
	trueInfos tInfos(aHit);
	
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	vector<G4double>      Edep = aHit->GetEdep();
	
	// Get scintillator volume x dimension (mm)
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.
	
	// Get Total Energy deposited
	double Etota = 0;
	double Ttota = 0;
	double latt  = 0;
	
	double att;
	
	double A    = ecc.attlen[sector-1][layer-1][0][strip-1];
	double B    = ecc.attlen[sector-1][layer-1][1][strip-1]*10.;
	double C    = ecc.attlen[sector-1][layer-1][2][strip-1];
	double G    = ecc.gain[sector-1][layer-1][strip-1];
	double a0   = ecc.timing[sector-1][layer-1][0][strip-1];
	double a2   = ecc.timing[sector-1][layer-1][2][strip-1];
	double veff = ecc.veff[sector-1][layer-1][strip-1]*10;
	
	for(unsigned int s=0; s<tInfos.nsteps; s++) {
		if(B>0) {
			double xlocal = Lpos[s].x();
			if(view==1) latt = pDx2 + xlocal;
			if(view==2) latt = pDx2 + xlocal;
			if(view==3) {
				if(layer > 3) {
					// for ec, it's a minus sign
					latt = pDx2-xlocal;
				} else {
					// for pcal, it's a plus sign
					latt = pDx2+xlocal;
				}
			}
			att   = A*exp(-latt/B)+C;
			Etota = Etota + Edep[s]*att;
			Ttota = Ttota + latt/veff;
		} else {
			Etota = Etota + Edep[s];
		}
	}
	
	// initialize ADC and TDC
	double ADC = 0;
	
	// simulate the adc value.
	if (Etota > 0) {
		double EC_npe = G4Poisson(Etota*pmtPEYld); //number of photoelectrons
		if (EC_npe>0) {
			double sigma  = sqrt(EC_npe)*ecc.pmtFactor;
			double EC_GeV = G4RandGauss::shoot(EC_npe,sigma)/1000./ecc.ADC_GeV_to_evio/G/pmtPEYld;
			if (EC_GeV>0) {
				ADC = EC_GeV;
				time_in_ns = (tInfos.time+Ttota/tInfos.nsteps) + a0 + a2/sqrt(ADC) + ecc.tdc_global_offset;
			}
		}
	}
	
	// Status flags
	if(accountForHardwareStatus) {
		switch (ecc.status[sector-1][layer-1][strip-1]) {
			case 0:
				break;
			case 1:
				ADC = 0;
				break;
			case 2:
				time_in_ns = 0;
				break;
			case 3:
				ADC = time_in_ns = 0;
				break;
				
			case 5:
				break;
				
			default:
				cout << " > Unknown EC status: " << ecc.status[sector-1][layer-1][strip-1] << " for sector " << sector << ",  layer " << layer << ", strip " << strip << endl;
		}		
	}
	
	// EVIO banks record time with offset determined by position of data in capture window.  On forward carriage this is currently
	// around 7.9 us.  This offset is omitted in the simulation.  Also EVIO TDC time is relative to the trigger time, which is not
	// simulated at present.
	
	int fadc_time = convert_to_precision(time_in_ns);
	int tdc = time_in_ns/a1;
	
	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = sector;
	dgtz["layer"]     = layer;
	dgtz["component"] = strip;
	dgtz["ADC_order"] = 0;
	dgtz["ADC_ADC"]   = ADC;
	dgtz["ADC_time"]  = fadc_time;
	dgtz["ADC_ped"]   = 0;
	dgtz["TDC_order"] = 2;
	dgtz["TDC_TDC"]   = tdc;
	
	// cout << "sector = " << sector << " layer = " << view << " strip = " << strip << " ADC = " << ADC << " TDC = " << TDC << endl;
	
	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}
	
	return dgtz;
}

vector<identifier>  ecal_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ecal_HitProcess :: electronicNoise()
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

map< string, vector <int> >  ecal_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// - charge: returns charge/time digitized information / step
// index 0: hit number
// index 1: step index
// index 2: charge
// index 3: time at electronics
// index 4: vector of identifiers - have to match the translation table
// index 5: hardware object: crate/slot/channel from translation table
map< int, vector <double> > ecal_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	vector<double> hitNumbers;
	vector<double> stepIndex;
	vector<double> chargeAtElectronics;
	vector<double> timeAtElectronics;
	vector<double> identifiers;
	vector<double> hardware;
	hitNumbers.push_back(hitn);
	
	// getting identifiers
	vector<identifier> identity = aHit->GetId();
	
	// get sector, stack (inner or outer), view (U, V, W), and strip.
	// The stack/view information is lost in the identifier but it can be recovered by looking at the volume description
	// For example:
	// grep "sector manual 1" ec__geometry_rga_fall2018.txt | awk -F\| '{print $2" "$18}' | grep "strip manual 1 "
	int sector = identity[0].id;
	int layer  = identity[1].id; // layer=1-3 (PCAL) 4-9 (ECAL). Layer = view for pcal, ecinner, ecouter
	int strip  = identity[2].id;
	
	// pcal
	int view   = layer;
	
	if (layer > 3 && layer < 7) {
		// ec inner (stack 1)
		view = layer - 3;
	} else if (layer > 7) {
		// ec inner (stack 2)
		view = layer - 6;
	}
	
	
	// Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
	// Different for EC and PCAL
	double pmtPEYld = 3.5 ;
	if (layer < 4) {
		// pcal
		pmtPEYld  = 11.5 ;
	}
	
	identifiers.push_back(sector);   // sector
	identifiers.push_back(layer);    // laylayer=1-3 (PCAL) 4-9 (ECAL)er
	identifiers.push_back(strip);    // component (pmt)
	identifiers.push_back(0);        // order
	
	// getting hardware
	Hardware thisHardware = ecc.TT.getHardware({sector, layer, strip, 0});
	hardware.push_back(thisHardware.getCrate());
	hardware.push_back(thisHardware.getSlot());
	hardware.push_back(thisHardware.getChannel());
	
	// Adding pedestal mean and sigma into the hardware as well
	// All of these variables start from 1, therefore -1 is subtracted, e.g. sector-1
	hardware.push_back(ecc.pedestal[sector - 1][layer - 1][view - 1]);
	hardware.push_back(ecc.pedestal_sigm[sector - 1][layer - 1][view - 1]);
	
	// getting charge and time
	trueInfos tInfos(aHit);
	
	// Get scintillator mother volume dimensions (mm)
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.
	
	vector<G4ThreeVector> pos  = aHit->GetPos();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	
	
	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> time = aHit->GetTime();
	
	double A  = ecc.attlen[sector-1][layer-1][0][strip-1];
	double B  = ecc.attlen[sector-1][layer-1][1][strip-1]*10.;
	double C  = ecc.attlen[sector-1][layer-1][2][strip-1];
	double G  = ecc.gain[sector-1][layer-1][strip-1];
	double veff  = ecc.veff[sector-1][layer-1][strip-1]*10;
	
	for(unsigned int s=0; s<tInfos.nsteps; s++) {
		if(B>0) {
			double xlocal = Lpos[s].x();
			double latt = 0;
			
			if(view==1) latt = pDx2+xlocal;
			if(view==2) latt = pDx2+xlocal;
			if(view==3) {
				if(layer > 3) {
					// for ecal, it's a minus sign
					latt = pDx2-xlocal;
				} else {
					// for ecal, it's a plus sign
					latt = pDx2+xlocal;
				}
			}
			double att   = A*exp(-latt/B)+C;
			
			double stepE = Edep[s]*att;
			double stepTime = time[s] + latt/veff;
			
			// cout<<"time[s] = "<<time[s]<<endl;
			// cout<<"att time  = "<<latt/ecc.veff<<endl;
			
			if (stepE > 0) {
				double EC_npe = G4Poisson(stepE*pmtPEYld); //number of photoelectrons
				if (EC_npe>0) {
					double sigma  = sqrt(EC_npe)*ecc.pmtFactor;
					double EC_GeV = G4RandGauss::shoot(EC_npe, sigma)/1000./ecc.ADC_GeV_to_evio/G/pmtPEYld;
					if (EC_GeV>0) {
						stepIndex.push_back(s);
						chargeAtElectronics.push_back(EC_GeV);
						timeAtElectronics.push_back(stepTime);
					}
				}
			}
		}
	}
	
	CT[0] = hitNumbers;
	CT[1] = stepIndex;
	CT[2] = chargeAtElectronics;
	CT[3] = timeAtElectronics;
	CT[4] = identifiers;
	CT[5] = hardware;
	
	
	// cout << " SIGNAL done with n steps: " << CT[3].size() << " and identifier " <<  CT[4].size() << endl;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double ecal_HitProcess :: voltage(double charge, double time, double forTime)
{
	//	return 0.0;
	// return DGauss(forTime, ecc.vpar, charge, time);
	return PulseShape(forTime, ecc.vpar, charge, time);
}

void ecal_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(ecc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ecc = initializeECConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		ecc.runNo = runno;
	}
}

// this static function will be loaded first thing by the executable
ecConstants ecal_HitProcess::ecc = initializeECConstants(-1);

