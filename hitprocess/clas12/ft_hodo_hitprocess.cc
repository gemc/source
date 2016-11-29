// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "ft_hodo_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static ftHodoConstants initializeFTHODOConstants(int runno)
{
	// all these constants should be read from CCDB
	ftHodoConstants fthc;


	// database
	fthc.runNo = runno;
    
    // do not initialize at the beginning, only after the end of the first event,
    // with the proper run number coming from options or run table
    if(runno == -1) return fthc;

	fthc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		fthc.connection = (string) getenv("CCDB_CONNECTION");
	else
		fthc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	fthc.variation  = "default";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(fthc.connection));

    
    int isector,ilayer,icomponent;

    vector<vector<double> > data;
    
    cout<<"FT-Hodo:Getting status"<<endl;
    sprintf(fthc.database,"/calibration/ft/fthodo/status:%d",fthc.runNo);
    data.clear(); calib->GetCalib(data,fthc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        isector    = data[row][0];
        ilayer     = data[row][1];
        icomponent = data[row][2];
        fthc.status[isector-1][ilayer-1].push_back(data[row][3]);
    }
    
    cout<<"FT-Hodo:Getting noise"<<endl;
    sprintf(fthc.database,"/calibration/ft/fthodo/noise:%d",fthc.runNo);
    data.clear(); calib->GetCalib(data,fthc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        isector    = data[row][0];
        ilayer     = data[row][1];
        icomponent = data[row][2];
        fthc.pedestal[isector-1][ilayer-1].push_back(data[row][3]);
        fthc.pedestal_rms[isector-1][ilayer-1].push_back(data[row][4]);
        fthc.gain_pc[isector-1][ilayer-1].push_back(data[row][5]);
        fthc.gain_mv[isector-1][ilayer-1].push_back(data[row][6]);
        fthc.npe_threshold[isector-1][ilayer-1].push_back(data[row][7]);
    }
    
    cout<<"FT-Hodo:Getting charge_to_energy"<<endl;
    sprintf(fthc.database,"/calibration/ft/fthodo/charge_to_energy:%d",fthc.runNo);
    data.clear(); calib->GetCalib(data,fthc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        isector    = data[row][0];
        ilayer     = data[row][1];
        icomponent = data[row][2];
        fthc.mips_charge[isector-1][ilayer-1].push_back(data[row][3]);
        fthc.mips_energy[isector-1][ilayer-1].push_back(data[row][4]);
    }
    
    cout<<"FT-Hodo:Getting time_offsets"<<endl;
    sprintf(fthc.database,"/calibration/ft/fthodo/time_offsets:%d",fthc.runNo);
    data.clear(); calib->GetCalib(data,fthc.database);
    for(unsigned row = 0; row < data.size(); row++)
    {
        isector    = data[row][0];
        ilayer     = data[row][1];
        icomponent = data[row][2];
        fthc.time_offset[isector-1][ilayer-1].push_back(data[row][3]);
        fthc.time_rms[isector-1][ilayer-1].push_back(data[row][4]);
    }

    // fadc parameters
    fthc.ns_per_sample = 4*ns;
    fthc.time_to_tdc   = 100/fthc.ns_per_sample;// conversion factor from time(ns) to TDC channels)
    fthc.tdc_max       = 8191;               // TDC range
    fthc.fadc_input_impedence = 50;
    fthc.fadc_LSB      = 0.4884;
    
    
    // setting voltage signal parameters
    fthc.vpar[0] = 20;  // delay, ns
    fthc.vpar[1] = 10;  // rise time, ns
    fthc.vpar[2] = 30;  // fall time, ns
    fthc.vpar[3] = 1;   // amplifier
    
    return fthc;
}

map<string, double> ft_hodo_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	trueInfos tInfos(aHit);
	
	// use Crystal ID to define IDX and IDY
	int isector    = identity[0].id;
	int ilayer     = identity[1].id;
    int icomponent = identity[2].id;
    
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = 8191;

	if(tInfos.eTot>0)
	{
        // adding shift and spread on time
        double time=tInfos.time+fthc.time_offset[isector-1][ilayer-1][icomponent-1]+G4RandGauss::shoot(0., fthc.time_rms[isector-1][ilayer-1][icomponent-1]);
        TDC=int(time*fthc.time_to_tdc);
        if(TDC>fthc.tdc_max) TDC=(int)fthc.tdc_max;

        // calculate charge and amplitude
        double charge    = tInfos.eTot*fthc.mips_charge[isector-1][ilayer-1][icomponent-1]/fthc.mips_energy[isector-1][ilayer-1][icomponent-1];
        double npe_mean  = charge/fthc.gain_pc[isector-1][ilayer-1][icomponent-1];
        double npe       = G4Poisson(npe_mean);
        charge           = charge * npe/npe_mean;
//        double amplitude = charge*fthc.gain_mv[isector-1][ilayer-1][icomponent-1]/fthc.gain_pc[isector-1][ilayer-1][icomponent-1];
//        double fadc      = amplitude/fthc.fadc_LSB;
        ADC = (int) (charge*fthc.fadc_input_impedence/fthc.fadc_LSB/fthc.ns_per_sample);
        
	}
	  
    // Status flags
    switch (fthc.status[isector-1][ilayer-1][icomponent-1])
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
            cout << " > Unknown FTHODO status: " << fthc.status[isector-1][ilayer-1][icomponent-1] << " for sector, layer, component "
                                                                                                   << isector << ", "
                                                                                                   << ilayer  << ", "
                                                                                                   << icomponent << endl;
    }

    dgtz["hitn"]      = hitn;
	dgtz["sector"]    = isector;
	dgtz["layer"]     = ilayer;
    dgtz["component"] = icomponent;
	dgtz["adc"]       = ADC;
	dgtz["tdc"]       = TDC;
		
	return dgtz;
}

vector<identifier>  ft_hodo_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


void ft_hodo_HitProcess::initWithRunNumber(int runno)
{
	if(fthc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		fthc = initializeFTHODOConstants(runno);
		fthc.runNo = runno;
	}
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ft_hodo_HitProcess :: electronicNoise()
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



map< string, vector <int> >  ft_hodo_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}


// - charge: returns charge/time digitized information / step
map< int, vector <double> > ft_hodo_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double ft_hodo_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

// this static function will be loaded first thing by the executable
ftHodoConstants ft_hodo_HitProcess::fthc = initializeFTHODOConstants(-1);








