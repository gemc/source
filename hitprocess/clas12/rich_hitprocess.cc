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

static richConstants initializeRICHConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
        richConstants richc;
	
	if(runno == -1) return richc;

        string timestamp = "";
        if(digiSnapshotTime != "no") {
                timestamp = ":"+digiSnapshotTime;
        }
	
	// database
	richc.runNo = runno;
	
	//	richc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != nullptr)
	  richc.connection = (string) getenv("CCDB_CONNECTION");
	else
	  richc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
        richc.variation  = "main";
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(richc.connection));

	vector<vector<double>> data;
	data.clear();
	// Eventually: not limited to module1

	// read timewalk correction
	snprintf(richc.database, sizeof(richc.database), "/calibration/rich/module1/time_walk:%d:%s%s", richc.runNo, digiVariation.c_str(), timestamp.c_str());
	calib->GetCalib(data,richc.database);
	
	for(int row = 0; row<data.size(); row++){	  
	  int ipmt = data[row][1];
	  richc.timewalkCorr_D0[ipmt-1] = data[row][3];
	  richc.timewalkCorr_m1[ipmt-1]	= data[row][4];
	  richc.timewalkCorr_m2[ipmt-1]	= data[row][5];
	  richc.timewalkCorr_T0[ipmt-1]	= data[row][6];	  
	}	
	//cout << "timewalk test: " << richc.timewalkCorr_D0[0] << endl;
	data.clear();
	  
        // read time offset
        snprintf(richc.database, sizeof(richc.database), "/calibration/rich/module1/time_offset:%d:%s%s", richc.runNo, digiVariation.c_str(), timestamp.c_str());
        calib->GetCalib(data,richc.database);

        for(int row = 0; row<data.size(); row++){
          int ipmt = data[row][1];
          richc.timeOffsetCorr[ipmt-1] = data[row][3];
        }

	// initialize rich pixel class
	richc.richPixel = new RichPixel();
	richc.richPixel->InitReadout(12700, 3e6, 1, 1);
	
	return richc;
}


// digitized info integrated over hit
// should this match what was in rich_sector4/rich__bank.txt?
// or should it roughly match RICH::tdc from data.json?

map<string, double> rich_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	//cout << "hit number " << hitn << endl;
        vector<identifier> identity = aHit->GetId();
        int idsector = identity[0].id;
        int idpmt = identity[1].id;
        int idpixel = 0;//identity[2].id;
	
	rejectHitConditions = false;
	writeHit = true;

	int pid  = aHit->GetPID();
	double stepzerotime = aHit->GetTime()[0];
	//cout << "step zero time: " << stepzerotime << endl;
	// setting all values generically just to see how they get printed out
	
	dgtz["hitn"]   = 0;//hitn;
	dgtz["sector"] = 0;//idsector;
	dgtz["layer"] = 0;//idpmt;
	dgtz["component"] = 0;//idpixel; // TODO: function for local position -> pixel #
	
	return dgtz;
}


// this routine needs to be modified
// no point drawing should be made here, but in MHit
// finding the PMT should be in a different function,
// with parameters coming from DB

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"

vector<identifier> rich_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
        // id[0]: sector
        // id[1]: pad
        // id[2]: pixel
        // id[i].id: number
        // time same for all
        // taken from dc_hitprocess: 
        G4StepPoint   *prestep   = aStep->GetPreStepPoint();
        G4StepPoint   *poststep  = aStep->GetPostStepPoint();
        G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction                          
	G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction                           
        ->GetTopTransform().TransformPoint(xyz);

	id[2].id = getPixelNumber(Lxyz);

	/*
	cout << "id size " << id.size() << endl;
	for(int i = 0; i <id.size(); i++){
	  cout << "detector identifier " << id[i].name << " " << id[i].time << " " << id[i].id << endl;  
	}
	*/
	//cout << "local pmt hit pos: " << Lxyz.x() << " " << Lxyz.y() << " " << Lxyz.z() << endl;
	//cout << "pixel number " << id[2].id << endl;
	id[id.size()-1].id_sharing = 1;
	return id;
}

// setting TDC information (need leading and trailing edge in rich reco)
map< string, vector <int> >  rich_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	vector<identifier> identity = aHit->GetId();
        int idsector = identity[0].id;
        int idpmt = identity[1].id;
        int idpixel = 0;//identity[2].id;

	int pid  = aHit->GetPID();
        double stepzerotime = aHit->GetTime()[0];

	// getting digitized timing information from RichPixel
	richc.richPixel->GenerateTDC(1, stepzerotime);
	// t2 and t1 now set

	vector<int> order, tdc;

	// leading edge
	order.push_back(2);
	tdc.push_back(richc.richPixel->get_T1());
	//trailing edge
	order.push_back(3);
	tdc.push_back(richc.richPixel->get_T2());	

	//cout << "t1 : " << richc.richPixel->get_T1() << endl;
	//cout <<	"t2 : "	<< richc.richPixel->get_T2() << endl;
	
	MH["TDC_TDC"]=tdc;
	MH["TDC_order"]=order;
        
        writeHit = true;	
	
	return MH;
}


void rich_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(richc.runNo != runno) {
		//		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		richc = initializeRICHConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		richc.runNo = runno;
	}
}

// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> rich_HitProcess :: electronicNoise()
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


// - charge: returns charge/time digitized information / step
map< int, vector <double> > rich_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The input are charge value, time
double rich_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 1.0;
}



// this static function will be loaded first thing by the executable
richConstants rich_HitProcess::richc = initializeRICHConstants(-1);

// PMT local position to pixel number
int rich_HitProcess::getPixelNumber(G4ThreeVector  Lxyz){
  // H8500:
  // 6.08mm for small
  // 6.26mm for large
  // H12700:
  // 6mm for small
  // 6.25mm for large
  // test: treating all as H12700
  // Pixel 1 is top left: -max x, +max y ?
  double edge_large = 6.25;
  double edge_small = 6.;
  
  double xloc = Lxyz.x() + 24.5; //mm
  double yloc = Lxyz.y() + 24.5;
  //cout << "x: " << xloc << " y: " << yloc << endl;
  int xpix = -1;
  int ypix = -1;
  if (xloc < edge_large){ xpix = 1; }
  else if (xloc > edge_large+6*edge_small){ xpix = 8; }
  else{ xpix = int((xloc-edge_large)/edge_small) + 1; }

  if (yloc < edge_large){ ypix = 1;}
  else if (yloc > edge_large+6*edge_small){ ypix = 8;}
  else{ ypix = int((yloc-edge_large)/edge_small) + 1; 
  }
  //cout << "xpix: " << xpix << " ypix: " << ypix << endl;
  return (int ((ypix-1)*8 + xpix));
}




/* ---------------------------------------------------*/
RichPixel::RichPixel(int t)
{

  if ( (t == 8500) || (t == 12700) || (t == 12701) ) {
    InitPmt(t);
  }
  else {
    printf("ERROR: pmt type %d not known\n", t);
    return;
  }


  InitMaroc();

  Clear();

  return;
}
/* ---------------------------------------------------*/
void RichPixel::InitPmt(int t, double g)
{

  PmtType = t;
  Gain = g;

  if (t == 8500) {
    nStages = 12;
    d1Ratio = 1;
  }
  else if ( (t == 12700) || (t == 12701) ) {
    nStages = 10;
    d1Ratio = 2;
  }
  else {
    printf("ERROR: pmt type %d not known\n", t);
    return;
  }


  GN = pow( (Gain/d1Ratio), 1./nStages);
  G1 = d1Ratio * GN;

  /* old calculation, wrong */
  if (t == 12701) {
    GN = pow(Gain, 1./(nStages+d1Ratio-1));
    G1 = pow(GN, d1Ratio);
  }


  return;
}
/* ---------------------------------------------------*/
void RichPixel::InitMaroc(double g, double th)
{

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
  for (int i=0; i<5; i++) {
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
void RichPixel::InitReadout(int pmttype, double pmtgain, double marocgain, double marocthr)
{
  InitPmt(pmttype, pmtgain);
  InitMaroc(marocgain, marocthr);

  return;
}
/* ---------------------------------------------------*/
void RichPixel::Clear()
{
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
int RichPixel::GenerateADC(int n0)
{

  GenerateNpe(n0);

  qadc = npe*Qe;
  ADC = Pedestal + (int)(DAC*qadc);

  return 0;
}
/* ---------------------------------------------------*/
bool RichPixel::GenerateTDC(int n0, double t0)
{

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
bool RichPixel::GenerateTDCfromADC(double qadc, double t0)
{

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
void RichPixel::GenerateNpe(int n0)
{
  int nEle = n0;
  for (int n=0; n<nStages; n++) {
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
void RichPixel::ChargeToTime()
{
  double qeff = qtdc / MarocThrF;
  double time = 0;
  if (qeff < q0_t) {
    for (int i=0; i<5; i++) {
      time = time + p_t[i] * pow(qeff, i);
    }
  }
  else {
    time = q_t + m_t * qeff;
  }

  true_t1 = time + TimeOffset + start_time;

  /* gaussian smearing of t1 */
  double dt = G4RandGauss::shoot(TimeOffset, TimeResol);

  t1 = time + dt + start_time;
  
  return;
}
/* -------------------------------------------------- */
void RichPixel::ChargeToDuration()
{
  //double qeff = qtdc / MarocThrF;
  double qeff = qtdc;
  double p0_eff = p0_d;
  if (MarocThrF > 1) {
    p0_eff = p0_d * (1. - (MarocThrF-1)*alphaD);
  }
  else if (MarocThrF < 1) {
    p0_eff = p0_d * (1. - (1./MarocThrF-1)*alphaD);
  }
  
  duration = p0_eff * (1. - exp( -sqrt(qeff/p1_d) ) );

  return;
}
/* -------------------------------------------------- */
void RichPixel::PrintPmt()
{

  printf("MAPMT H%d\n", PmtType);
  printf("  G=%f    nStages=%d   g1=%f   gN=%f\n", Gain, nStages, G1, GN);

  return;
}
/* -------------------------------------------------- */
void RichPixel::PrintMaroc()
{

  printf("MAROC setting\n");
  printf("  gain=%f   RelThreshold=%f   threshold=%f fC\n", MarocG, MarocThrF, MarocThrCharge);

  return;
}
