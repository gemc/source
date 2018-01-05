// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"
#include <math.h>

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "pcal_hitprocess.h"

static pcConstants initializePCConstants(int runno)
{
	pcConstants pcc;

	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return pcc;

	int isec,ilay;
//	int isec,ilay,istr;

	// database
	pcc.runNo      = runno;
	pcc.date       = "2015-11-29";

	if(getenv ("CCDB_CONNECTION") != NULL) {
	  pcc.connection = (string) getenv("CCDB_CONNECTION");
	} else {
	  pcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	}
	
	pcc.variation  = "default";

	pcc.TDC_time_to_evio    = 1.      ;  	                                   
	pcc.ADC_GeV_to_evio     = 1./10000; // MIP based calibration is nominally 10 channels/MeV
	pcc.veff                = 160.    ; // Effective velocity of scintillator light (mm/ns)
	pcc.pmtPEYld            = 11.5    ; // Number of p.e. divided by the energy deposited in MeV. See EC NIM paper table 1.
	pcc.pmtQE               = 0.27    ;
	pcc.pmtDynodeGain       = 4.0     ;
	pcc.pmtDynodeK          = 0.5     ; // K=0 (Poisson) K=1(exponential)	vector<vector<double> > data;
	//  Fluctuations in PMT gain distributed using Gaussian with
	//  sigma=1/SNR where SNR = sqrt[(1-QE+(k*del+1)/(del-1))/npe] del = dynode gain k=0-1
	//  Adapted from G-75 (pg. 169) and and G-111 (pg. 174) from RCA PMT Handbook.
	//  Factor k for dynode statistics can range from k=0 (Poisson) to k=1 (exponential).
	//  Note: GSIM sigma was incorrect (used 1/sigma for sigma).
	pcc.pmtFactor           = sqrt(1-pcc.pmtQE+(pcc.pmtDynodeK*pcc.pmtDynodeGain+1)/(pcc.pmtDynodeGain-1));

	vector<vector<double> > data;
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(pcc.connection));

	sprintf(pcc.database,"/calibration/ec/gain:%d",pcc.runNo);
	data.clear(); calib->GetCalib(data,pcc.database);
    
	for(unsigned row = 0; row < data.size(); row++)
	  {
	    isec = data[row][0]; ilay = data[row][1];
		  //istr = data[row][2];
	    pcc.gain[isec-1][ilay-1].push_back(data[row][3]);
	  }

	sprintf(pcc.database,"/calibration/ec/attenuation:%d",pcc.runNo);
	data.clear(); calib->GetCalib(data,pcc.database);

	for(unsigned row = 0; row < data.size(); row++)
	  {
	    isec = data[row][0]; ilay = data[row][1];
		  //istr = data[row][2];
	    pcc.attlen[isec-1][ilay-1][0].push_back(data[row][3]);
	    pcc.attlen[isec-1][ilay-1][1].push_back(data[row][5]);
	    pcc.attlen[isec-1][ilay-1][2].push_back(data[row][7]);
	  }
	
	sprintf(pcc.database,"/calibration/ec/timing:%d",pcc.runNo);
	data.clear(); calib->GetCalib(data,pcc.database);

	for(unsigned row = 0; row < data.size(); row++)
	  {
	    isec = data[row][0]; ilay = data[row][1];
		  //istr = data[row][2];
	    pcc.timing[isec-1][ilay-1][0].push_back(data[row][3]);
	    pcc.timing[isec-1][ilay-1][1].push_back(data[row][4]);
	    pcc.timing[isec-1][ilay-1][2].push_back(data[row][5]);
	    pcc.timing[isec-1][ilay-1][3].push_back(data[row][6]);
	    pcc.timing[isec-1][ilay-1][4].push_back(data[row][7]);
	  }
	
	// FOR now we will initialize pedestals and sigmas to a random value, in the future
	// they will be initialized from CCDB, when the DB will be ready
	const double const_ped_value = 101;
	const double const_ped_sigm_value = 2;
	// commands below fill all the elements of pcc.pedestal and pcc.pedestal_sigm with their values (const_ped_value, and const_ped_sigm_value respectively)
	std::fill(&pcc.pedestal[0][0][0], &pcc.pedestal[0][0][0] + sizeof(pcc.pedestal)/sizeof(pcc.pedestal[0][0][0]), const_ped_value);
	std::fill(&pcc.pedestal_sigm[0][0][0], &pcc.pedestal_sigm[0][0][0] + sizeof(pcc.pedestal_sigm)/sizeof(pcc.pedestal_sigm[0][0][0]), const_ped_sigm_value);
	
	// setting voltage signal parameters
	pcc.vpar[0] = 0;  // delay, ns
	pcc.vpar[1] = 2.8;  // rise time, ns
	pcc.vpar[2] = 20;  // fall time, ns
	pcc.vpar[3] = 1;   // amplifier


	// loading translation table
	pcc.TT = TranslationTable("pcTT");

	// loads translation table from CLAS12 Database:
	// Translation table for EC (ECAL+PCAL).
	// Crate sector assignments: ECAL/FADC=1,7,13,19,25,31 ECAL/TDC=2,8,14,20,26,32
	// PCAL/FADC=3,9,15,21,27,33 PCAL/TDC=4,10,16,22,28,34.
	// ORDER: 0=FADC 2=TDC.

	string database   = "/daq/tt/ec:1";


	data.clear(); calib->GetCalib(data, database);
	cout << "  > " << pcc.TT.getName() << " TT Data loaded from CCDB with " << data.size() << " columns." << endl;

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
		pcc.TT.addHardwareItem({sector, layer, pmt, order}, Hardware(crate, slot, channel));
	}
	cout << "  > Data loaded in translation table " << pcc.TT.getName() << endl;
	
	return pcc;
}


void pcal_HitProcess::initWithRunNumber(int runno)
{
	if(pcc.runNo != runno)
	{
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		pcc = initializePCConstants(runno);
		pcc.runNo = runno;
	}
}
map<string, double> pcal_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;

	vector<identifier> identity = aHit->GetId();

	// get sector, view (U, V, W), and strip.
	int sector = identity[0].id;
	int module = identity[1].id;
	int view   = identity[2].id;
	int strip  = identity[3].id;
	trueInfos tInfos(aHit);

	HCname = "PCAL Hit Process";

	// Get scintillator volume x dimension (mm)
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.

	//cout<<"pDx2="<<pDx2<<" pDy1="<<pDy1<<endl;


	// Get Total Energy deposited
	double Etota = 0;
	double Ttota = 0;
	double latt  = 0;

	vector<G4double>      Edep = aHit->GetEdep();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();

	double att;

	double A  = pcc.attlen[sector-1][view-1][0][strip-1];
	double B  = pcc.attlen[sector-1][view-1][1][strip-1]*10.;
	double C  = pcc.attlen[sector-1][view-1][2][strip-1];
        double G  = pcc.gain[sector-1][view-1][strip-1];
	double a0 = pcc.timing[sector-1][view-1][0][strip-1];
	double a1 = pcc.timing[sector-1][view-1][1][strip-1];
	double a2 = pcc.timing[sector-1][view-1][2][strip-1];
	double a3 = pcc.timing[sector-1][view-1][3][strip-1]/100;
	double a4 = pcc.timing[sector-1][view-1][4][strip-1]/1000;
	
	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		if(B>0)
		{
			double xlocal = Lpos[s].x();
			if(view==1) latt=pDx2+xlocal;
			if(view==2) latt=pDx2+xlocal;
			if(view==3) latt=pDx2+xlocal;
			att   = A*exp(-latt/B)+C;
			Etota = Etota + Edep[s]*att;
			Ttota = Ttota + latt/pcc.veff + a3*latt*latt + a4*latt*latt*latt;
		}
		else
		{
			Etota = Etota + Edep[s];
		}
	}

	// initialize ADC and TDC
	double ADC = 0;
	double TDC = 0;

	// simulate the adc value.
	if (Etota > 0) {
		double PC_npe = G4Poisson(Etota*pcc.pmtPEYld); //number of photoelectrons
		if (PC_npe>0) {
			double sigma = pcc.pmtFactor/sqrt(PC_npe);
			double PC_GeV = G4RandGauss::shoot(PC_npe,sigma)/1000./pcc.ADC_GeV_to_evio/G/pcc.pmtPEYld;
			if (PC_GeV>0) {
				ADC = PC_GeV;
				TDC = (tInfos.time+Ttota/tInfos.nsteps)*pcc.TDC_time_to_evio + a0 + a2/sqrt(ADC);
			}
		}
	}

	// EVIO banks record time with offset determined by position of data in capture window.  On forward carriage this is currently
	// around 7.9 us.  This offset is omitted in the simulation.  Also EVIO TDC time is relative to the trigger time, which is not
	// simulated at present.

	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["module"] = module;
	dgtz["view"]   = view;
	dgtz["strip"]  = strip;
	dgtz["ADC"]    = ADC;
	dgtz["TDC"]    = TDC/a1;

	//cout << "sector = " << sector << " layer = " << module << " view = " << view << " strip = " << strip << " PL_ADC = " << ADC << " TDC = " << TDC << " Edep = " << Etot << endl;

	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	bool rejectHitConditions = false;
	if(rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
}


vector<identifier>  pcal_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	//int sector             = yid[0].id;
	//int layer              = yid[1].id;
	//int view               = yid[2].id; // get the view: U->1, V->2, W->3
	//int strip	       = yid[3].id;
	//return yid;
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> pcal_HitProcess :: electronicNoise()
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


map< string, vector <int> >  pcal_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;

	return MH;
}



// - charge: returns charge/time digitized information / step
map< int, vector <double> > pcal_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	vector<double> hitNumbers;
	vector<double> stepIndex;
	vector<double> chargeAtElectronics;
	vector<double> timeAtElectronics;
	vector<double> identifiers;
	vector<double> hardware;
	hitNumbers.push_back(hitn);


	vector<identifier> identity = aHit->GetId();

	// get sector, stack (inner or outer), view (U, V, W), and strip.
	int sector = identity[0].id;
	int stack  = identity[1].id;
	int view   = identity[2].id;
	int strip  = identity[3].id;
	int layer  = (stack-1)*3+view; // layer=1-3 (PCAL) 4-9 (ECAL)

	identifiers.push_back(sector);
	identifiers.push_back(layer);
	identifiers.push_back(strip);    // component (pmt)
	identifiers.push_back(0);        // order

	// getting hardware
	Hardware thisHardware = pcc.TT.getHardware({sector, layer, strip, 0});
	hardware.push_back(thisHardware.getCrate());
	hardware.push_back(thisHardware.getSlot());
	hardware.push_back(thisHardware.getChannel());


	// Adding pedestal mean and sigma into the hardware as well
	// All of these variables start from 1, therefore -1 is sbutracted, e.g. sector-1
	hardware.push_back(pcc.pedestal[sector - 1][layer - 1][view - 1]);
	hardware.push_back(pcc.pedestal_sigm[sector - 1][layer - 1][view - 1]);


	trueInfos tInfos(aHit);
	
	// Get scintillator mother volume dimensions (mm)
	//double pDy1 = aHit->GetDetector().dimensions[3];  ///< G4Trap Semilength.
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.
	//double BA   = sqrt(4*pow(pDy1,2) + pow(pDx2,2)) ;

	vector<G4ThreeVector> pos  = aHit->GetPos();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();


	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> time = aHit->GetTime();

	double A  = pcc.attlen[sector-1][layer-1][0][strip-1];
	double B  = pcc.attlen[sector-1][layer-1][1][strip-1]*10.;
	double C  = pcc.attlen[sector-1][layer-1][2][strip-1];
	double G  = pcc.gain[sector-1][layer-1][strip-1];

	for(unsigned int s=0; s<tInfos.nsteps; s++) {
		if(B>0) {
			double xlocal = Lpos[s].x();
			//double ylocal = Lpos[s].y();
			double latt = 0;

			//if(view==1) latt = xlocal+(pDx2/(2.*pDy1))*(ylocal+pDy1);
			//if(view==2) latt = BA*(pDy1-ylocal)/2./pDy1;
			//if(view==3) latt = BA*(ylocal+pDy1-xlocal*2*pDy1/pDx2)/4/pDy1;
		      	if(view==1) latt = pDx2+xlocal;
	       		if(view==2) latt = pDx2+xlocal;
			if(view==3) latt = pDx2-xlocal;
			double att   = A*exp(-latt/B)+C;

			double stepE = Edep[s]*att;
			double stepTime = time[s] + latt/pcc.veff;

			// cout<<"time[s] = "<<time[s]<<endl;
			// cout<<"att time  = "<<latt/pcc.veff<<endl;

			if (stepE > 0) {
				double PC_npe = G4Poisson(stepE*pcc.pmtPEYld); //number of photoelectrons
				if (PC_npe>0) {
					double sigma  = pcc.pmtFactor/sqrt(PC_npe);
					double PC_GeV = G4RandGauss::shoot(PC_npe, sigma)/1000./pcc.ADC_GeV_to_evio/G/pcc.pmtPEYld;
					if (PC_GeV>0) {
						stepIndex.push_back(s);
						chargeAtElectronics.push_back(PC_GeV);
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
double pcal_HitProcess :: voltage(double charge, double time, double forTime)
{
	//	return 0.0;
//	return DGauss(forTime, pcc.vpar, charge, time);
	return PulseShape(forTime, pcc.vpar, charge, time);
}


// this static function will be loaded first thing by the executable
pcConstants pcal_HitProcess::pcc = initializePCConstants(-1);







