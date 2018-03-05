// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// gemc headers
#include "ftof_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

static ftofConstants initializeFTOFConstants(int runno) {
    ftofConstants ftc;

    // do not initialize at the beginning, only after the end of the first event,
    // with the proper run number coming from options or run table
    if (runno == -1) return ftc;

    ftc.runNo = runno;
    ftc.date = "2015-11-29";
    if (getenv("CCDB_CONNECTION") != NULL)
        ftc.connection = (string) getenv("CCDB_CONNECTION");
    else
        ftc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
    ftc.variation = "default";

    ftc.npaddles[0] = 23;
    ftc.npaddles[1] = 62;
    ftc.npaddles[2] = 5;

    ftc.thick[0] = 5.0;
    ftc.thick[1] = 6.0;
    ftc.thick[2] = 5.0;

    ftc.dEdxMIP = 1.956; // muons in polyvinyltoluene
    ftc.pmtPEYld = 500;
    //	ftc.tdcLSB        = 42.5532;// counts per ns (23.5 ps LSB)

    cout << "FTOF:Setting time resolution" << endl;
    for (int p = 0; p < 3; p++) {
        for (int c = 1; c < ftc.npaddles[p] + 1; c++) {
            if (p == 0) ftc.tres[p].push_back(1e-3 * (c * 5.45 + 74.55)); //ps to ns
            if (p == 1) ftc.tres[p].push_back(1e-3 * (c * 0.90 + 29.10)); //ps to ns
            if (p == 2) ftc.tres[p].push_back(1e-3 * (c * 5.00 + 145.0)); //ps to ns
        }
    }

    ftc.dEMIP[0] = ftc.thick[0] * ftc.dEdxMIP;
    ftc.dEMIP[1] = ftc.thick[1] * ftc.dEdxMIP;
    ftc.dEMIP[2] = ftc.thick[2] * ftc.dEdxMIP;

    int isec, ilay, istr;

    vector<vector<double> > data;

    auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(ftc.connection));
    cout << "Connecting to " << ftc.connection << "/calibration/ftof" << endl;

    cout << "FTOF:Getting attenuation" << endl;
    sprintf(ftc.database, "/calibration/ftof/attenuation:%d", ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.attlen[isec - 1][ilay - 1][0].push_back(data[row][3]);
        ftc.attlen[isec - 1][ilay - 1][1].push_back(data[row][4]);
    }

    cout << "FTOF:Getting effective_velocity" << endl;
    sprintf(ftc.database, "/calibration/ftof/effective_velocity:%d", ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.veff[isec - 1][ilay - 1][0].push_back(data[row][3]);
        ftc.veff[isec - 1][ilay - 1][1].push_back(data[row][4]);
    }

    cout << "FTOF:Getting status" << endl;
    sprintf(ftc.database, "/calibration/ftof/status:%d", ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.status[isec - 1][ilay - 1][0].push_back(data[row][3]);
        ftc.status[isec - 1][ilay - 1][1].push_back(data[row][4]);
    }

    cout << "FTOF:Getting gain_balance" << endl;
    sprintf(ftc.database, "/calibration/ftof/gain_balance:%d", ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.countsForMIP[isec - 1][ilay - 1][0].push_back(data[row][3]);
        ftc.countsForMIP[isec - 1][ilay - 1][1].push_back(data[row][4]);
    }

    cout << "FTOF:Getting time_walk" << endl;
    sprintf(ftc.database, "/calibration/ftof/time_walk:%d", ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.twlk[isec - 1][ilay - 1][0].push_back(data[row][3]);
        ftc.twlk[isec - 1][ilay - 1][1].push_back(data[row][4]);
        ftc.twlk[isec - 1][ilay - 1][2].push_back(data[row][5]);
        ftc.twlk[isec - 1][ilay - 1][3].push_back(data[row][6]);
        ftc.twlk[isec - 1][ilay - 1][4].push_back(data[row][7]);
        ftc.twlk[isec - 1][ilay - 1][5].push_back(data[row][8]);
    }

    cout << "FTOF:Getting time_offset" << endl;
    
    sprintf(ftc.database,"/calibration/ftof/time_offsets:%d",ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.toff_LR[isec - 1][ilay - 1].push_back(data[row][3]);
        ftc.toff_RFpad[isec-1][ilay-1].push_back(data[row][4]);
        ftc.toff_P2P[isec-1][ilay-1].push_back(data[row][5]);
    
    }

    cout << "FTOF:Getting tdc_conv" << endl;
    sprintf(ftc.database, "/calibration/ftof/tdc_conv:%d", ftc.runNo);
    data.clear();
    calib->GetCalib(data, ftc.database);
    for (unsigned row = 0; row < data.size(); row++) {
        isec = data[row][0];
        ilay = data[row][1];
        istr = data[row][2];
        ftc.tdcconv[isec - 1][ilay - 1][0].push_back(data[row][3]);
        ftc.tdcconv[isec - 1][ilay - 1][1].push_back(data[row][4]);
    }

    // setting voltage signal parameters
    ftc.vpar[0] = 50; // delay, ns
    ftc.vpar[1] = 10; // rise time, ns
    ftc.vpar[2] = 20; // fall time, ns
    ftc.vpar[3] = 1; // amplifier


    // FOR now we will initialize pedestals and sigmas to a random value, in the future
    // they will be initialized from DB 
    const double const_ped_value = 101;
    const double const_ped_sigm_value = 2;
    // commands below fill all the elements of ctc.pedestal and ctc.pedestal_sigm with their values (const_ped_value, and const_ped_sigm_value respectively)
    std::fill(&ftc.pedestal[0][0][0][0], &ftc.pedestal[0][0][0][0] + sizeof (ftc.pedestal) / sizeof (ftc.pedestal[0][0][0][0]), const_ped_value);
    std::fill(&ftc.pedestal_sigm[0][0][0][0], &ftc.pedestal_sigm[0][0][0][0] + sizeof (ftc.pedestal_sigm) / sizeof (ftc.pedestal_sigm[0][0][0][0]), const_ped_sigm_value);

    string database = "/daq/tt/ftof:1";

    data.clear();
    calib->GetCalib(data, database);
    cout << "  > " << ftc.TT.getName() << " TT Data loaded from CCDB with " << data.size() << " columns." << endl;

    // filling translation table
    for (unsigned row = 0; row < data.size(); row++) {
        int crate = data[row][0];
        int slot = data[row][1];
        int channel = data[row][2];

        int sector = data[row][3];
        int panel = data[row][4];
        int paddle = data[row][5];
        int pmt = data[row][6];

        // order is important as we could have duplicate entries w/o it
        ftc.TT.addHardwareItem({sector, panel, paddle, pmt}, Hardware(crate, slot, channel));

    }
    cout << "  > Data loaded in translation table " << ftc.TT.getName() << endl;

    return ftc;
}

void ftof_HitProcess::initWithRunNumber(int runno) {
    if (ftc.runNo != runno) {
        cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
        ftc = initializeFTOFConstants(runno);
        ftc.runNo = runno;
    }
}

map<string, double> ftof_HitProcess::integrateDgt(MHit* aHit, int hitn) {
    map<string, double> dgtz;
	if(aHit->isBackgroundHit == 1) return dgtz;

    vector<identifier> identity = aHit->GetId();

    int sector = identity[0].id;
    int panel = identity[1].id;
    int paddle = identity[2].id;
    int pmt = identity[3].id; // 0=> Left PMT, 1=> Right PMT

    trueInfos tInfos(aHit);

    // Get the paddle half-length 
    double length = aHit->GetDetector().dimensions[0];

    // Distances from left, right
    //	double dLeft  = length + tInfos.lx;
    //	double dRight = length - tInfos.lx;

    double d = length + (1 - 2 * pmt) * tInfos.lx;

    // attenuation length
    //	double attlenL = ftc.attlen[sector-1][panel-1][0][paddle-1];
    //	double attlenR = ftc.attlen[sector-1][panel-1][1][paddle-1];

    double attlen = ftc.attlen[sector - 1][panel - 1][pmt][paddle - 1];
    double attlen_otherside = ftc.attlen[sector - 1][panel - 1][1 - pmt].at(paddle - 1);
    // attenuation factor
    //	double attLeft  = exp(-dLeft/cm/attlenL);
    //	double attRight = exp(-dRight/cm/attlenR);

    double att = exp(-d / cm / attlen);

    // Gain factors to simulate FTOF PMT gain matching algorithm.
    // Each L,R PMT pair has HV adjusted so geometeric mean sqrt(L*R)
    // is independent of counter length, which compensates for
    // the factor exp(-L/2/attlen) where L=full length of bar.
    //	double gainLeft  = sqrt(attLeft*attRight);
    //	double gainRight = gainLeft;

    double gain = sqrt(exp(-d / cm / attlen) * exp(-(2 * length - d) / cm / attlen_otherside));

    // Attenuated light at PMT
    //	double eneL = tInfos.eTot*attLeft;
    //	double eneR = tInfos.eTot*attRight;

    double ene = tInfos.eTot*att;

    // TDC conversion factors
    double tdcconv = ftc.tdcconv[sector - 1][panel - 1][pmt][paddle - 1];

    // giving geantinos some energies
    if (aHit->GetPID() == 0) {
        double gmomentum = aHit->GetMom().mag() / GeV;
        //		eneL = gmomentum*attLeft;
        //		eneR = gmomentum*attRight;

        ene = gmomentum*att;

    }

    double adc = 0;
    double adcu = 0;
    double tdc = 0;
    double tdcu = 0;

    // Fluctuate the light measured by the PMT with
    // Poisson distribution for emitted photoelectrons
    // Treat L and R separately, in case nphe=0

    if (ene > 0) {
        adcu = ene * ftc.countsForMIP[sector - 1][panel - 1][pmt][paddle - 1] / ftc.dEMIP[panel - 1] / gain;
    }



    double nphe = G4Poisson(ene * ftc.pmtPEYld);
    ene = nphe / ftc.pmtPEYld;

    if (ene > 0) {
        adc = ene * ftc.countsForMIP[sector - 1][panel - 1][pmt][paddle - 1] / ftc.dEMIP[panel - 1] / gain;
        double A = ftc.twlk[sector - 1][panel - 1][3 * pmt + 0][paddle - 1];
        double B = ftc.twlk[sector - 1][panel - 1][3 * pmt + 1][paddle - 1];
        //double            C = ftc.twlk[sector-1][panel-1][2][paddle-1];
        double timeWalk = A / pow(adc, B);
//        double tU = tInfos.time + d / ftc.veff[sector - 1][panel - 1][pmt][paddle - 1] / cm +
//                ftc.toff_LR[sector - 1][panel - 1][paddle - 1] / 2. - ftc.toff_P2P[sector - 1][panel - 1][paddle - 1] + timeWalk;
        
        double tU = tInfos.time + d/ftc.veff[sector-1][panel-1][pmt][paddle-1]/cm + ftc.toff_LR[sector-1][panel-1][paddle-1]/2. 
                                                                                                - ftc.toff_RFpad[sector-1][panel-1][paddle-1] 
                                                                                                - ftc.toff_P2P[sector-1][panel-1][paddle-1] 
                                                                                                + timeWalk; 
        
        
        double t = G4RandGauss::shoot(tU, sqrt(2) * ftc.tres[panel - 1][paddle - 1]);
        tdcu = tU / tdcconv;
        tdc = t / tdcconv;
    }


    // Status flags
    switch (ftc.status[sector - 1][panel - 1][pmt][paddle - 1]) {
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
            cout << " > Unknown FTOF status: " << ftc.status[sector - 1][panel - 1][0][paddle - 1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " left " << endl;
    }


    //	cout << " > FTOF status: " << ftc.status[sector-1][panel-1][0][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " left: " << adcl << endl;
    //	cout << " > FTOF status: " << ftc.status[sector-1][panel-1][1][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " right:  " << adcr << endl;


    dgtz["hitn"] = hitn;
    dgtz["sector"] = sector;
    dgtz["layer"] = panel;
    dgtz["paddle"] = paddle;
    dgtz["side"] = (int) pmt;
    dgtz["ADC"] = (int) adc;
    dgtz["TDC"] = (int) tdc;
    dgtz["ADCu"] = (int) adcu;
    dgtz["TDCu"] = (int) tdcu;

    // decide if write an hit or not
    writeHit = true;
    // define conditions to reject hit
    bool rejectHitConditions = false;
    if (rejectHitConditions) {
        writeHit = false;
    }

    return dgtz;
}

vector<identifier> ftof_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector) {

    id[id.size() - 1].id_sharing = 1;

    vector<identifier> yid = id;
    yid[0].id_sharing = 1; // sector
    yid[1].id_sharing = 1; // panel
    yid[2].id_sharing = 1; // paddle
    yid[3].id_sharing = 1; // side, left or right

    if (yid[3].id != 0) {
        cout << "*****WARNING***** in ftof_HitProcess :: processID, identifier PTT of the original hit should be 0 " << endl;
        cout << "yid[3].id = " << yid[3].id << endl;
    }

    // Now we want to have similar identifiers, but the only difference be id PMT to be 1, instead of 0    
    identifier this_id = yid[0];
    yid.push_back(this_id);
    this_id = yid[1];
    yid.push_back(this_id);
    this_id = yid[2];
    yid.push_back(this_id);
    this_id = yid[3];
    this_id.id = 1;
    yid.push_back(this_id);

    return yid;
}

// - electronicNoise: returns a vector of hits generated / by electronics.

vector<MHit*> ftof_HitProcess::electronicNoise() {
    vector<MHit*> noiseHits;

    // first, identify the cells that would have electronic noise
    // then instantiate hit with energy E, time T, identifier IDF:
    //
    // MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

    // push to noiseHits collection:
    // noiseHits.push_back(thisNoiseHit)

    return noiseHits;
}

map< string, vector <int> > ftof_HitProcess::multiDgt(MHit* aHit, int hitn) {
    map< string, vector <int> > MH;

    return MH;
}

// - charge: returns charge/time digitized information / step

map< int, vector <double> > ftof_HitProcess::chargeTime(MHit* aHit, int hitn) {
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

    int sector = identity[0].id;
    int panel = identity[1].id;
    int paddle = identity[2].id;
    int pmt = identity[3].id; // 0=> Left PMT, 1=> Right PMT

    identifiers.push_back(sector); // sector
    identifiers.push_back(panel); // panel, 1a, 1b, 2a
    identifiers.push_back(paddle); // paddle number
    identifiers.push_back(pmt); // the pmt side: 0=> Left, 1=>Right

    // getting hardware
    Hardware thisHardware = ftc.TT.getHardware({sector, panel, paddle, pmt});
    hardware.push_back(thisHardware.getCrate());
    hardware.push_back(thisHardware.getSlot());
    hardware.push_back(thisHardware.getChannel());

    // Adding pedestal mean and sigma into the hardware as well
    // All of these variables start from 1, therefore -1 is subtracted, e.g. sector-1
    hardware.push_back(ftc.pedestal[sector - 1][panel - 1][paddle - 1][pmt]);
    hardware.push_back(ftc.pedestal_sigm[sector - 1][panel - 1 ][paddle - 1][pmt]);

    // attenuation length
    double attlen = ftc.attlen[sector - 1][panel - 1][pmt][paddle - 1];
    double attlen_otherside = ftc.attlen[sector - 1][panel - 1][1 - pmt].at(paddle - 1);

    trueInfos tInfos(aHit);

    // Get the paddle half-length 
    double length = aHit->GetDetector().dimensions[0];

    // Vector of positions of the hit in each step
    vector<G4ThreeVector> Lpos = aHit->GetLPos();

    // Vector of Edep and time of the hit in each step
    vector<G4double> Edep = aHit->GetEdep();
    vector<G4double> time = aHit->GetTime();

    for (unsigned int s = 0; s < tInfos.nsteps; s++) {
        // Distances from left, right
        //	double dLeft  = length + tInfos.lx;
        //	double dRight = length - tInfos.lx;
        
        double d = length + (1 - 2 * pmt) * Lpos[s].x();
        //double d = length + (1 - 2 * pmt) * tInfos.lx;

        // attenuation factor
        //	double attLeft  = exp(-dLeft/cm/attlenL);
        //	double attRight = exp(-dRight/cm/attlenR);

        double att = exp(-d / cm / attlen);

        // Gain factors to simulate FTOF PMT gain matching algorithm.
        // Each L,R PMT pair has HV adjusted so geometeric mean sqrt(L*R)
        // is independent of counter length, which compensates for
        // the factor exp(-L/2/attlen) where L=full length of bar.
        //	double gainLeft  = sqrt(attLeft*attRight);
        //	double gainRight = gainLeft;

        double gain = sqrt(exp(-d / cm /attlen ) * exp(-(2 * length - d) / cm / attlen_otherside));

        // Attenuated light at PMT
        //	double eneL = tInfos.eTot*attLeft;
        //	double eneR = tInfos.eTot*attRight;

        double ene = Edep[s] * att;

        // giving geantinos some energies
        if (aHit->GetPID() == 0) {
            double gmomentum = aHit->GetMom().mag() / GeV;
            //		eneL = gmomentum*attLeft;
            //		eneR = gmomentum*attRight;

            ene = gmomentum*att;
        }

        double adc = 0;
        double adcu = 0;

		// Fluctuate the light measured by the PMT with
        // Poisson distribution for emitted photoelectrons
        // Treat L and R separately, in case nphe=0

        if (ene > 0) {
            adcu = ene * ftc.countsForMIP[sector - 1][panel - 1][pmt][paddle - 1] / ftc.dEMIP[panel - 1] / gain;
        }

        double nphe = G4Poisson(ene * ftc.pmtPEYld);
        ene = nphe / ftc.pmtPEYld;

        if (ene > 0) {
            adc = ene * ftc.countsForMIP[sector - 1][panel - 1][pmt][paddle - 1] / ftc.dEMIP[panel - 1] / gain;
            double A = ftc.twlk[sector - 1][panel - 1][3 * pmt + 0][paddle - 1];
            double B = ftc.twlk[sector - 1][panel - 1][3 * pmt + 1][paddle - 1];
            //double            C = ftc.twlk[sector-1][panel-1][2][paddle-1];
            double timeWalk = A / pow(adc, B);
            
            double stepTimeU = time[s] + d/ftc.veff[sector-1][panel-1][pmt][paddle-1]/cm + ftc.toff_LR[sector-1][panel-1][paddle-1]/2. 
                                                                                                - ftc.toff_RFpad[sector-1][panel-1][paddle-1] 
                                                                                                - ftc.toff_P2P[sector-1][panel-1][paddle-1] 
                                                                                                + timeWalk;
            
            double stepTime = G4RandGauss::shoot(stepTimeU, sqrt(2) * ftc.tres[panel - 1][paddle - 1]);

            stepIndex.push_back(s); // Since it is going to be only one hit, i.e. only one step
            chargeAtElectronics.push_back(adc);
            timeAtElectronics.push_back(stepTime);
        }

    }

    //	// Status flags
    //	switch (ftc.status[sector-1][panel-1][pmt][paddle-1])
    //	{
    //		case 0:
    //			break;
    //		case 1:
    //			adc = 0;
    //			break;
    //		case 2:
    //			tdc = 0;
    //			break;
    //		case 3:
    //			adc = tdc = 0;
    //			break;
    //
    //		case 5:
    //			break;
    //			
    //		default:
    //			cout << " > Unknown FTOF status: " << ftc.status[sector-1][panel-1][0][paddle-1] << " for sector " << sector << ",  panel " << panel << ", paddle " << paddle << " left " << endl;
    //	}


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

double ftof_HitProcess::voltage(double charge, double time, double forTime) {
    //	return 0.0;
    return PulseShape(forTime, ftc.vpar, charge, time);
}

// this static function will be loaded first thing by the executable
ftofConstants ftof_HitProcess::ftc = initializeFTOFConstants(-1);




