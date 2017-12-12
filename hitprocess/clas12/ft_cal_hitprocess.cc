// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "ft_cal_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;



static ftCalConstants initializeFTCALConstants(int runno)
{
	// all these constants should be read from CCDB
	ftCalConstants ftcc;

    // do not initialize at the beginning, only after the end of the first event,
    // with the proper run number coming from options or run table
    if(runno == -1) return ftcc;

	// database
	ftcc.runNo = runno;
	ftcc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		ftcc.connection = (string) getenv("CCDB_CONNECTION");
	else
		ftcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	ftcc.variation  = "default";

    int icomponent;
    
    vector<vector<double> > data;

    auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ftcc.connection));
    cout<<"Connecting to "<<ftcc.connection<<"/calibration/ft/ftcal"<<endl;

    cout<<"FT-Cal:Getting status"<<endl;
    sprintf(ftcc.database,"/calibration/ft/ftcal/status:%d",ftcc.runNo);
    data.clear(); calib->GetCalib(data,ftcc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        icomponent   = data[row][2];
        ftcc.status[icomponent] = data[row][3];
    }
    
    cout<<"FT-Cal:Getting noise"<<endl;
    sprintf(ftcc.database,"/calibration/ft/ftcal/noise:%d",ftcc.runNo);
    data.clear(); calib->GetCalib(data,ftcc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        icomponent   = data[row][2];
        ftcc.pedestal[icomponent] = data[row][3];         ftcc.pedestal[icomponent] = 101.;   // When DB will be filled, I should remove this
        ftcc.pedestal_rms[icomponent] = data[row][4];     ftcc.pedestal_rms[icomponent] = 2.; // When DB will be filled, I should remove this
        ftcc.noise[icomponent] = data[row][5];
        ftcc.noise_rms[icomponent] = data[row][6];
        ftcc.threshold[icomponent] = data[row][7];
    }

    cout<<"FT-Cal:Getting charge_to_energy"<<endl;
    sprintf(ftcc.database,"/calibration/ft/ftcal/charge_to_energy:%d",ftcc.runNo);
    data.clear(); calib->GetCalib(data,ftcc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        icomponent   = data[row][2];
        ftcc.mips_charge[icomponent] = data[row][3];
        ftcc.mips_energy[icomponent] = data[row][4];
        ftcc.fadc_to_charge[icomponent] = data[row][5];
        ftcc.preamp_gain[icomponent] = data[row][6];
        ftcc.apd_gain[icomponent] = data[row][7];
    }

    cout<<"FT-Cal:Getting time_offsets"<<endl;
    sprintf(ftcc.database,"/calibration/ft/ftcal/time_offsets:%d",ftcc.runNo);
    data.clear(); calib->GetCalib(data,ftcc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        icomponent   = data[row][2];
        ftcc.time_offset[icomponent] = data[row][3];
        ftcc.time_rms[icomponent] = data[row][4];
    }

    
    cout<<"FT-Cal: Getting Translation table"<<endl;
    
    	// loading translation table
	ftcc.TT = TranslationTable("ftcalTT");

	// loads translation table from CLAS12 Database:
	// Translation table for ft cal.
	// Sector is always 1, layer is 0, component is the channel number, and order is always 0.

	string database   = "/daq/tt/ftcal:1";


	data.clear(); calib->GetCalib(data, database);
	cout << "  > " << ftcc.TT.getName() << " TT Data loaded from CCDB with " << data.size() << " columns." << endl;

	// filling translation table
	for(unsigned row = 0; row < data.size(); row++)
	{
		int crate   = data[row][0];
		int slot    = data[row][1];
		int channel = data[row][2];

		int sector  = data[row][3];
		int layer   = data[row][4];
		int crystal = data[row][5];
		int order   = data[row][6];

		// order is important as we could have duplicate entries w/o it
		ftcc.TT.addHardwareItem({sector, layer, crystal, order}, Hardware(crate, slot, channel));
	}
	cout << "  > Data loaded in translation table " << ftcc.TT.getName() << endl;

    
    
    // fadc parameters
    ftcc.ns_per_sample = 4*ns;
    ftcc.time_to_tdc   = 100/ftcc.ns_per_sample;// conversion factor from time(ns) to TDC channels)
    ftcc.tdc_max       = 8191;               // TDC range

    // preamp parameters
    ftcc.preamp_input_noise = 5500;     // preamplifier input noise in number of electrons
    ftcc.apd_noise          = 0.0033;   // relative noise based on a Voltage and Temperature stability of 10 mV (3.9%/V) and 0.1 C (3.3%/C)
    
    // crystal paramters
    ftcc.light_speed = 15*cm/ns;
    
    // setting voltage signal parameters
    ftcc.vpar[0] = 20;  // delay, ns
    ftcc.vpar[1] = 10;  // rise time, ns
    ftcc.vpar[2] = 30;  // fall time, ns
    ftcc.vpar[3] = 1;   // amplifier
    
    
    return ftcc;
}

map<string, double> ft_cal_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;

	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);

	
	// R.De Vita (November 2016)
	
/* old digitization
	// relevant parameter for digitization (in the future should be read from database)
	double adc_charge_tochannel=20;    // conversion factor from charge(pC) to ADC channels
	double PbWO4_light_yield =240/MeV; // Lead Tungsten Light Yield (APD have similar QE for
	// fast component, lambda=420nm-ly=120ph/MeV, and slow component,
	// lambda=560nm-ly=20ph/MeV, taking fast component only)
	// double PbWO4_light_yield =672/MeV; // LY at -25 deg=2.8 x LY at +18 deg
	double APD_qe    = 0.70;           // APD Quantum Efficiency (Hamamatsu S8664-55)
	double APD_size  = 100*mm*mm;       // APD size ( 10 mm x 10 mm)
	double APD_noise = 0.0033;          // relative noise based on a Voltage and Temperature stability of 10 mV (3.9%/V) and 0.1 C (3.3%/C)
 
    double mips_energy = 15.3*MeV;
    double mips_charge = 6.005;
    double ns_per_sample = 4*ns;
    double fadc_input_impedence = 50;
    double fadc_LSB = 0.4884;
    double fadc_to_charge = fadc_LSB*ns_per_sample/fadc_input_impedence;
	
    double time_to_tdc=100/ns_per_sample;// conversion factor from time(ns) to TDC channels)
    double tdc_max=8191;               // TDC range
    double time_res=0.2;               // time resolution
    double light_speed =15;

    double AMP_input_noise = 5500;     // preamplifier input noise in number of electrons
    double AMP_gain        = 700;      // preamplifier gain
    double APD_gain  = 150;            // based on FT note
 */

    // Get the crystal length: in the FT crystal are BOXes and the half-length is the 3rd element
    double length = 2 * aHit->GetDetector().dimensions[2];
	// Get the crystal width (rear face): in the FT crystal are BOXes and the half-length is the 2th element
    //	double width  = 2 * aHit->GetDetector().dimensions[1];
	
	// use Crystal ID to define IDX and IDY
	int IDX = identity[0].id;
	int IDY = identity[1].id;
    
        int iCrystal = (IDY-1)*22+IDX-1;
    
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = 8191;
	
	
	
	if(tInfos.eTot>0)
	{
		/*
		 // commented out to use average time instead of minimum time in TDC calculation
		 for(int s=0; s<nsteps; s++)
		 {
     double dRight = length/2 - Lpos[s].z();              // distance along z between the hit position and the end of the crystal
     double timeR  = times[s] + dRight/cm/light_speed;    // arrival time of the signal at the end of the crystal (speed of light in the crystal=15 cm/ns)
     if(Edep[s]>1*MeV) Tmin=min(Tmin,timeR);              // signal time is set to first hit time with energy above 1 MeV
		 }
		 TDC=int(Tmin*tdc_time_to_channel);
		 if(TDC>tdc_max) TDC=(int)tdc_max;
		 */
		double dRight = length/2 - tInfos.lz;                 // distance along z between the hit position and the end of the crystal
		double timeR  = tInfos.time + dRight/ftcc.light_speed;  // arrival time of the signal at the end of the crystal (speed of light in the crystal=15 cm/ns)
		// adding shift and spread on time
		timeR=timeR+ftcc.time_offset[iCrystal]+G4RandGauss::shoot(0., ftcc.time_rms[iCrystal]);
		
		TDC=int(timeR*ftcc.time_to_tdc);
		if(TDC>ftcc.tdc_max) TDC=(int)ftcc.tdc_max;
		
		// calculate number of photoelectrons detected by the APD considering the light yield, the q.e., and the size of the sensor
        //old!!		double npe=G4Poisson(tInfos.eTot*PbWO4_light_yield*0.5*APD_qe*APD_size/width/width);
        // for PMT, an addition factor of 0.5625 is needed to reproduce the 13.5 photoelectrons with a 20% QE
        //   double npe=G4Poisson(Etot*PbWO4_light_yield*0.5*0.5625*APD_qe*APD_size/width/width);
        double charge   = tInfos.eTot*ftcc.mips_charge[iCrystal]/ftcc.mips_energy[iCrystal];

        // add spread due to photoelectron statistics
        double npe_mean = charge/1.6E-7/ftcc.preamp_gain[iCrystal]/ftcc.apd_gain[iCrystal];
        double npe      = G4Poisson(npe_mean);
        
		// calculating APD output charge (in number of electrons) and adding noise
		double nel=npe*ftcc.apd_gain[iCrystal];
		nel=nel*G4RandGauss::shoot(1.,ftcc.apd_noise);
		if(nel<0) nel=0;
		// adding preamplifier input noise
		nel=nel+ftcc.preamp_input_noise*G4RandGauss::shoot(0.,1.);
		if(nel<0) nel=0;
        
		// converting to charge (in picoCoulomb)
        //old!		double crg=nel*AMP_gain*1.6e-7;
        //        double crg = charge * nel/(npe_mean*ftcc.apd_gain[iCrystal]);
        ADC = (int) (charge/ftcc.fadc_to_charge[iCrystal]);
        //old! ADC= (int) (crg*adc_charge_tochannel);
   	
	}
    
    // Status flags
    switch (ftcc.status[iCrystal])
    {
        case 0:
            break;
        case 1:
            break;
        case 3:
            ADC = TDC = 0;
            break;
            
        case 5:
            break;
            
        default:
            cout << " > Unknown FTCAL status: " << ftcc.status[iCrystal] << " for component " << iCrystal << endl;
    }

	
	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = 1;
    dgtz["layer"]     = 1;
	dgtz["component"] = iCrystal;
	dgtz["adc"]       = ADC;
	dgtz["tdc"]       = TDC;
	
	// decide if write an hit or not
	writeHit = true;
	// define conditions to reject hit
	bool rejectHitConditions = false;
	if(rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
}

vector<identifier>  ft_cal_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

void ft_cal_HitProcess::initWithRunNumber(int runno)
{
	if(ftcc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ftcc = initializeFTCALConstants(runno);
		ftcc.runNo = runno;
	}
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ft_cal_HitProcess :: electronicNoise()
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



map< string, vector <int> >  ft_cal_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}




// - charge: returns charge/time digitized information / step
map< int, vector <double> > ft_cal_HitProcess :: chargeTime(MHit* aHit, int hitn) {
    map< int, vector <double> > CT;

    vector<double> hitNumbers;
    vector<double> stepIndex;
    vector<double> chargeAtElectronics;
    vector<double> timeAtElectronics;
    vector<double> identifiers;
    vector<double> hardware;
    hitNumbers.push_back(hitn);

    // getting identifiers
    vector<identifier> identity = aHit->GetId();

    int sector = 1; // Always 1
    int layer = 0; // Always 0;
    int iX = identity[0].id;
    int iY = identity[1].id;
    int crystal = (iY - 1) * 22 + iX - 1;
    int order = 0; // Always 0

    identifiers.push_back(sector); // sector Always 1
    identifiers.push_back(layer); // layer Always 0
    identifiers.push_back(crystal); // component This is the crystal index
    identifiers.push_back(order); // Always 0

    // getting hardware
    Hardware thisHardware = ftcc.TT.getHardware({sector, layer, crystal, 0});
    hardware.push_back(thisHardware.getCrate());
    hardware.push_back(thisHardware.getSlot());
    hardware.push_back(thisHardware.getChannel());

    // Adding pedestal mean and sigma into the hardware as well
    // All of these variables start from 1, therefore -1 is subtracted, e.g. sector-1
    hardware.push_back(ftcc.pedestal[crystal]);
    hardware.push_back(ftcc.pedestal_rms[crystal]);


    trueInfos tInfos(aHit);

    // Get the crystal length: in the FT crystal are BOXes and the half-length is the 3rd element
    double length = 2 * aHit->GetDetector().dimensions[2];

    vector<G4ThreeVector> Lpos = aHit->GetLPos();

    vector<G4double> Edep = aHit->GetEdep();
    vector<G4double> time = aHit->GetTime();
    
    for (unsigned int s = 0; s < tInfos.nsteps; s++) {

        double dRight = length / 2 - Lpos[s].z(); // distance along z between the hit position and the end of the crystal
        double stepTime = time[s] + dRight / ftcc.light_speed; // arrival time of the signal at the end of the crystal (speed of light in the crystal=15 cm/ns)
        // adding shift and spread on time
        stepTime = stepTime + ftcc.time_offset[crystal] + G4RandGauss::shoot(0., ftcc.time_rms[crystal]);

        //cout<<"stepTime = "<<stepTime<<endl;
//        if( stepTime > 400 ){
//            //cout<<"==================== STEPTIME > 400 ns"<<endl;
//            cout<<"Step = "<<s<<"   stepTime = "<<stepTime<<"     Energy = "<<Edep[s]<<endl;
//        }
            
        
        // calculate number of photoelectrons detected by the APD considering the light yield, the q.e., and the size of the sensor
        //old!!		double npe=G4Poisson(tInfos.eTot*PbWO4_light_yield*0.5*APD_qe*APD_size/width/width);
        // for PMT, an addition factor of 0.5625 is needed to reproduce the 13.5 photoelectrons with a 20% QE
        //   double npe=G4Poisson(Etot*PbWO4_light_yield*0.5*0.5625*APD_qe*APD_size/width/width);
        double stepCharge = Edep[s] * ftcc.mips_charge[crystal] / ftcc.mips_energy[crystal];

        // add spread due to photoelectron statistics
        double npe_mean = stepCharge / 1.6E-7 / ftcc.preamp_gain[crystal] / ftcc.apd_gain[crystal];
        double npe = G4Poisson(npe_mean);

        // calculating APD output charge (in number of electrons) and adding noise
        double nel = npe * ftcc.apd_gain[crystal];
        nel = nel * G4RandGauss::shoot(1., ftcc.apd_noise);
        if (nel < 0) nel = 0;
        // adding preamplifier input noise
        nel = nel + ftcc.preamp_input_noise * G4RandGauss::shoot(0., 1.);
        if (nel < 0) nel = 0;

        // converting to charge (in picoCoulomb)
        //old!		double crg=nel*AMP_gain*1.6e-7;
        //        double crg = charge * nel/(npe_mean*ftcc.apd_gain[iCrystal]);
        
        // It is better to have it double at tis moment, it will be converted to in the MeventAction
        double ADC = (stepCharge / ftcc.fadc_to_charge[crystal]);
        //old! ADC= (int) (crg*adc_charge_tochannel);

        stepIndex.push_back(s);
        chargeAtElectronics.push_back(ADC);
        timeAtElectronics.push_back(stepTime);



    }
        
    // === Testing =======
//    double tot_adc = 0;
//    for( int ii = 0; ii < chargeAtElectronics.size(); ii++ ){
//        tot_adc = tot_adc + chargeAtElectronics.at(ii);
//        cout<<"cur_adc and time = "<<chargeAtElectronics.at(ii)<<"   "<<timeAtElectronics.at(ii);
//    }
//    cout<<"\n Tot adc = "<<tot_adc<<endl;
    
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
double ft_cal_HitProcess :: voltage(double charge, double time, double forTime)
{
	return PulseShape(forTime, ftcc.vpar, charge, time);
}


// this static function will be loaded first thing by the executable
ftCalConstants ft_cal_HitProcess::ftcc = initializeFTCALConstants(-1);






