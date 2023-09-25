// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>

using namespace ccdb;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"

using namespace CLHEP;

// gemc headers
#include "rich_hitprocess.h"

static richConstants
initializeRICHConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false) {
    // TODO: with TDC simulation class from Marco M., time calibration information maybe not necessary
    richConstants richc;
    if (runno == -1) return richc;

    string timestamp = "";
    if (digiSnapshotTime != "no") {
        timestamp = ":" + digiSnapshotTime;
    }

    // database
    richc.runNo = runno;

    //	richc.date       = "2016-03-15";
    if (getenv("CCDB_CONNECTION") != nullptr)
        richc.connection = (string) getenv("CCDB_CONNECTION");
    else
        richc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

    richc.variation = "main";
    unique_ptr <Calibration> calib(CalibrationGenerator::CreateCalibration(richc.connection));

    vector <vector<double>> data;
    data.clear();
    // Eventually: not limited to module1

    // read timewalk correction
    snprintf(richc.database, sizeof(richc.database), "/calibration/rich/module1/time_walk:%d:%s%s", richc.runNo, digiVariation.c_str(), timestamp.c_str());
    calib->GetCalib(data, richc.database);

    for (unsigned int row = 0; row < data.size(); row++) {
        int ipmt = data[row][1];
        richc.timewalkCorr_D0[ipmt - 1] = data[row][3];
        richc.timewalkCorr_m1[ipmt - 1] = data[row][4];
        richc.timewalkCorr_m2[ipmt - 1] = data[row][5];
        richc.timewalkCorr_T0[ipmt - 1] = data[row][6];
    }

    data.clear();

    // read time offset
    snprintf(richc.database, sizeof(richc.database), "/calibration/rich/module1/time_offset:%d:%s%s", richc.runNo, digiVariation.c_str(), timestamp.c_str());
    calib->GetCalib(data, richc.database);

    for (unsigned int row = 0; row < data.size(); row++) {
        int ipmt = data[row][1];
        richc.timeOffsetCorr[ipmt - 1] = data[row][3];
    }

    return richc;
}


// digitized info integrated over hit
// changed to match data.json definition of RICH::tdc 

map<string, double> rich_HitProcess::integrateDgt(MHit *aHit, int hitn) {
    map<string, double> dgtz;

    trueInfos tInfos(aHit);

    vector <identifier> identity = aHit->GetId();
    vector<double> time = aHit->GetTime();
    int idsector = identity[0].id;

    // tdc bank expects tile number
    int idpmt = identity[1].id;
    int tile = richc.pmtToTile[idpmt - 1];

    // tdc bank: readout channel number
    // pixel, order, tdc already set in processID
    int idpixel = identity[2].userInfos[2];//identity[2].id;
    int marocChannel = richc.anodeToMaroc[idpixel - 1];
    int tileChannel = marocChannel + (richc.pmtToTilePosition[idpmt - 1] - 1) * 64;

    int order = identity[2].userInfos[0];
    int tdc = identity[2].userInfos[1];

    writeHit = true;
    rejectHitConditions = false;

    double energy = aHit->GetEs()[0] / electronvolt;

    double qeff = 0;
    for (int i = 0; i < richc.nQEbins; i++) {
        if (energy < richc.Ene[i] && energy > richc.Ene[i + 1]) {
            qeff = richc.QE[i];
            break;
        }
    }

    // applying quantum efficiency from thrown random value set in integrateDgt
    if (identity[2].userInfos[3] > qeff) {
        writeHit = false;
    }

    // only time offset effects, duration calculated in intDgt,
    // so adding hit time to tdc here
    dgtz["hitn"] = hitn;
    dgtz["sector"] = idsector;
    dgtz["layer"] = tile;
    dgtz["component"] = tileChannel;
    dgtz["TDC_TDC"] = tdc + int(time[0]);
    dgtz["TDC_order"] = order;
    return dgtz;
}

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"

vector <identifier> rich_HitProcess::processID(vector <identifier> id, G4Step *aStep, detector Detector) {
    vector <identifier> yid = id;
    // id[0]: sector
    // id[1]: pmt
    // id[2]: pixel

    G4StepPoint *prestep = aStep->GetPreStepPoint();
    G4StepPoint *poststep = aStep->GetPostStepPoint();
    G4ThreeVector xyz = poststep->GetPosition();                                        ///< Global Coordinates of interaction
    G4ThreeVector Lxyz = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
            ->GetTopTransform().TransformPoint(xyz);

    int pixel = getPixelNumber(Lxyz);

    G4ThreeVector pixelCenterLocal = getPixelCenter(pixel);
    G4ThreeVector pixelCenterGlobal = prestep->GetTouchableHandle()->GetHistory()->GetTopTransform().Inverse().TransformPoint(pixelCenterLocal);


    int pmt = yid[1].id;
    RichPixel richPixel(richc.pmtType[pmt - 1]);

    richPixel.Clear();

    int t1 = -1;
    int t2 = -1;

    if (richPixel.GenerateTDC(1, 0)) {
        // generating TDC
        t1 = richPixel.get_T1();
        t2 = richPixel.get_T2();
    }

    for (int i = 0; i < 3; i++) {
        identifier idtemp;
        idtemp.name = id[i].name;
        idtemp.rule = id[i].rule;
        idtemp.id = id[i].id;
        idtemp.TimeWindow = id[i].TimeWindow;
        idtemp.TrackId = id[i].TrackId;

        yid.push_back(idtemp);
    }

    double QEthrow = G4UniformRand();
    yid[2].id = pixel;

    yid[2].userInfos.clear();
    yid[2].userInfos.push_back(double(1)); // TDC_order
    yid[2].userInfos.push_back(double(t1)); // TDC_tdc
    yid[2].userInfos.push_back(double(pixel)); // pixel
    yid[2].userInfos.push_back(QEthrow); // thrown random value for quantum eff.

    yid[5].id = pixel;
    yid[5].userInfos.clear();
    yid[5].userInfos.push_back(double(0));
    yid[5].userInfos.push_back(double(t2));
    yid[5].userInfos.push_back(double(pixel));
    yid[5].userInfos.push_back(QEthrow);

    yid[2].id_sharing = .5;
    yid[5].id_sharing = .5;

    return yid;
}

// setting TDC information (need leading and trailing edge in rich reco)
map <string, vector<int>> rich_HitProcess::multiDgt(MHit *aHit, int hitn) {
    map <string, vector<int>> MH;
    return MH;
}


void rich_HitProcess::initWithRunNumber(int runno) {
    string digiVariation = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
    string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;

    if (richc.runNo != runno) {
        //		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
        richc = initializeRICHConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
        richc.runNo = runno;
    }
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit *> rich_HitProcess::electronicNoise() {
    vector < MHit * > noiseHits;

    // first, identify the cells that would have electronic noise
    // then instantiate hit with energy E, time T, identifier IDF:
    //
    // MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

    // push to noiseHits collection:
    // noiseHits.push_back(thisNoiseHit)

    return noiseHits;
}


// - charge: returns charge/time digitized information / step
map<int, vector<double> > rich_HitProcess::chargeTime(MHit *aHit, int hitn) {
    map<int, vector<double> > CT;

    return CT;
}

// - voltage: returns a voltage value for a given time. The input are charge value, time
double rich_HitProcess::voltage(double charge, double time, double forTime) {
    return 1.0;
}


// this static function will be loaded first thing by the executable
richConstants rich_HitProcess::richc = initializeRICHConstants(-1);


// PMT local position to pixel number
int rich_HitProcess::getPixelNumber(G4ThreeVector Lxyz) {
    // H8500:
    // 6.08mm for small
    // 6.26mm for large
    // H12700:
    // 6mm for small
    // 6.25mm for large
    // test: treating all as H12700
    // Pixel 1 is top left: -max x, +max y ?
   // double edge_large = 6.25;
    double edge_small = 6.;

    double xloc = Lxyz.x(); //mm
    double yloc = Lxyz.y();

    int nx = int(abs(xloc / edge_small));
    int ny = int(abs(yloc / edge_small));
    int direcx = (xloc / abs(xloc));
    int direcy = (yloc / abs(yloc));

    if (nx > 3) nx = 3;
    if (nx < -3) nx = -3;
    if (ny > 3) ny = 3;
    if (ny < -3) ny = -3;

    int xpix = 0;
    int ypix = 0;
    if (nx > 0 && direcx == 1) {
        xpix = -1 * nx + 4;
    }
    if (nx == 0 && direcx == 1) {
        xpix = 4;
    }
    if (nx == 0 && direcx == -1) {
        xpix = 5;
    }
    if (nx > 0 && direcx == -1) {
        xpix = (nx + 1) + 4;
    }


    if (ny > 0 && direcy == 1) {
        ypix = -1 * ny + 4;
    }
    if (ny == 0 && direcy == 1) {
        ypix = 4;
    }
    if (ny == 0 && direcy == -1) {
        ypix = 5;
    }
    if (ny > 0 && direcy == -1) {
        ypix = (ny + 1) + 4;
    }

    return (int((ypix - 1) * 8 + xpix));
}

G4ThreeVector rich_HitProcess::getPixelCenter(int pixel) {
    // center is (0,0)
    // 1 should be -x, -y
    int xpix = int((pixel - 1) % 8) + 1;
    int ypix = int((pixel - 1) / 8) + 1;
    double edge_large = 6.25;
    double edge_small = 6.;

    double xpos = 0;
    double ypos = 0;

    if (xpix == 1) {
        xpos += (3 * edge_small + 0.5 * edge_large);
    }
    if (xpix == 8) {
        xpos -= (3 * edge_small + 0.5 * edge_large);
    }

    if (ypix == 1) {
        ypos += (3 * edge_small + 0.5 * edge_large);
    }
    if (ypix == 8) {
        ypos -= (3 * edge_small + 0.5 * edge_large);
    }

    if (xpix > 1 && xpix <= 4) {
        xpos += -1 * ((xpix - 4) - 0.5) * edge_small;
    }
    if (xpix >= 5 && xpix < 8) {
        xpos += -1 * (xpix - 5 + 0.5) * edge_small;
    }

    if (ypix > 1 && ypix <= 4) {
        ypos += -1 * ((ypix - 4) - 0.5) * edge_small;
    }
    if (ypix >= 5 && ypix < 8) {
        ypos += -1 * (ypix - 5 + 0.5) * edge_small;
    }

    return G4ThreeVector(xpos, ypos, -0.05);


}

/* ---------------------------------------------------*/
RichPixel::RichPixel(int t) {

    if ((t == 8500) || (t == 12700) || (t == 12701)) {
        InitPmt(t);
    } else {
        printf("ERROR: pmt type %d not known\n", t);
        return;
    }


    InitMaroc();

    Clear();

    return;
}

/* ---------------------------------------------------*/
void RichPixel::InitPmt(int t, double g) {

    PmtType = t;
    Gain = g;

    if (t == 8500) {
        nStages = 12;
        d1Ratio = 1;
    } else if ((t == 12700) || (t == 12701)) {
        nStages = 10;
        d1Ratio = 2;
    } else {
        printf("ERROR: pmt type %d not known\n", t);
        return;
    }


    GN = pow((Gain / d1Ratio), 1. / nStages);
    G1 = d1Ratio * GN;

    /* old calculation, wrong */
    if (t == 12701) {
        GN = pow(Gain, 1. / (nStages + d1Ratio - 1));
        G1 = pow(GN, d1Ratio);
    }


    return;
}

/* ---------------------------------------------------*/
void RichPixel::InitMaroc(double g, double th) {

    MarocG = g;
    MarocThrF = th;

    ThresholdScale = 5;

    MarocThrCharge = ThresholdScale * (MarocThrF - 1) * StdTrhesholdDAC * DAC;



    /* charge to time conversion, saturated region */
    q0_t = 230;
    /* Values from Ctest data analysis
      p_t[0] = 118.3;
    p_t[1] = -0.3134;
    p_t[2] = 0.002787;
    p_t[3] = -1.382e-5;
    p_t[4] = 2.531e-8;
    */

    p_t[0] = 7.32 + 4;
    p_t[1] = -0.08027;
    p_t[2] = 0.0001407;
    p_t[3] = 0;
    p_t[4] = 0;


    /* charge to time conversion, linear region */
    //m_t = -0.00307; //value from CTest data analysis
    m_t = -0.00207;
    q_t = -m_t * q0_t;
    for (int i = 0; i < 5; i++) {
        q_t = q_t + p_t[i] * pow(q0_t, i);
    }


    /* from charge to duration */

    //p0_d = 80; //from CTest data
    //p1_d = 170.1;  //from CTest data

    p0_d = 68.19;
    p1_d = 60.1;

    /* parameter to shift doen the duration per 1 unit of MarocThrF */
    alphaD = 0.1;


    return;
}

/* ---------------------------------------------------*/
void RichPixel::InitReadout(int pmttype, double pmtgain, double marocgain, double marocthr) {
    InitPmt(pmttype, pmtgain);
    InitMaroc(marocgain, marocthr);

    return;
}

/* ---------------------------------------------------*/
void RichPixel::Clear() {
    npe = 0;

    qadc = 0;
    ADC = 0;

    qtdc = 0;
    start_time = 0;
    true_t1 = 0;
    t1 = 0;
    t2 = 0;
    duration = 0;


    return;
}

/* ---------------------------------------------------*/
int RichPixel::GenerateADC(int n0) {

    GenerateNpe(n0);

    qadc = npe * Qe;
    ADC = Pedestal + (int) (DAC * qadc);

    return 0;
}

/* ---------------------------------------------------*/
bool RichPixel::GenerateTDC(int n0, double t0) {
    GenerateNpe(n0);
    qtdc = npe * Qe * MarocG;
    if (qtdc > MarocMaxQ) qtdc = MarocMaxQ;

    if (qtdc < MarocThrCharge) return false;

    start_time = t0;
    ChargeToTime();
    ChargeToDuration();
    t2 = t1 + duration;

    return true;
}

/* ---------------------------------------------------*/
bool RichPixel::GenerateTDCfromADC(double qadc, double t0) {

    qtdc = qadc * MarocG;
    if (qtdc > MarocMaxQ) qtdc = MarocMaxQ;

    if (qtdc < MarocThrCharge) return false;

    start_time = t0;
    ChargeToTime();
    ChargeToDuration();
    t2 = t1 + duration;


    return true;
}

/* -------------------------------- */
void RichPixel::GenerateNpe(int n0) {
    int nEle = n0;
    for (int n = 0; n < nStages; n++) {
        double g = GN;
        if (n == 0) g = G1;

        int nIn = nEle;
        double nAve = g * nIn;
        nEle = G4Poisson(nAve);
        if (nEle == 0) nEle = n0;
    }
    npe = nEle;


    return;
}

/* -------------------------------------------------- */
void RichPixel::ChargeToTime() {
    double qeff = qtdc / MarocThrF;
    double time = 0;
    if (qeff < q0_t) {
        for (int i = 0; i < 5; i++) {
            time = time + p_t[i] * pow(qeff, i);
        }
    } else {
        time = q_t + m_t * qeff;
    }

    true_t1 = time + TimeOffset + start_time;

    /* gaussian smearing of t1 */
    double dt = G4RandGauss::shoot(TimeOffset, TimeResol);

    t1 = time + dt + start_time;

    return;
}

/* -------------------------------------------------- */
void RichPixel::ChargeToDuration() {
    //double qeff = qtdc / MarocThrF;
    double qeff = qtdc;
    double p0_eff = p0_d;
    if (MarocThrF > 1) {
        p0_eff = p0_d * (1. - (MarocThrF - 1) * alphaD);
    } else if (MarocThrF < 1) {
        p0_eff = p0_d * (1. - (1. / MarocThrF - 1) * alphaD);
    }

    duration = p0_eff * (1. - exp(-sqrt(qeff / p1_d)));

    return;
}

/* -------------------------------------------------- */
void RichPixel::PrintPmt() {

    printf("MAPMT H%d\n", PmtType);
    printf("  G=%f    nStages=%d   g1=%f   gN=%f\n", Gain, nStages, G1, GN);

    return;
}

/* -------------------------------------------------- */
void RichPixel::PrintMaroc() {

    printf("MAROC setting\n");
    printf("  gain=%f   RelThreshold=%f   threshold=%f fC\n", MarocG, MarocThrF, MarocThrCharge);

    return;
}
