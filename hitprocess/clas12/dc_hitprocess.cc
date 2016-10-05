// gemc headers
#include "dc_hitprocess.h"

// clhep
#include "CLHEP/Random/RandGaussT.h"
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// geant4
#include "Randomize.hh"


static dcConstants initializeDCConstants(int runno)
{
	// all these constants should be read from CCDB
	dcConstants dcc;

	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return dcc;
	
	dcc.NWIRES = 113;
	dcc.dcThreshold  = 50;  // eV

	// database
	dcc.runNo = runno;
	dcc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		dcc.connection = (string) getenv("CCDB_CONNECTION");
	else
		dcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	dcc.variation  = "main";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(dcc.connection));

	
	// reading efficiency parameters
	string database   = "/calibration/dc/signal_generation/intrinsic_inefficiency";
	vector<vector<double> > data;
	calib->GetCalib(data, database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		int sl = data[row][0] - 1;
		dcc.P1[sl]     = data[row][1];
		dcc.P2[sl]     = data[row][2];
		dcc.P3[sl]     = data[row][3];
		dcc.P4[sl]     = data[row][4];
		dcc.iScale[sl] = data[row][5];
	}
	
	// reading smearing parameters
	database   = "/calibration/dc/signal_generation/dc_resolution";
	data.clear();
	calib->GetCalib(data, database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		int sec = data[row][0] - 1;
		int sl  = data[row][1] - 1;
		dcc.smearP1[sec][sl]    = data[row][2];
		dcc.smearP2[sec][sl]    = data[row][3];
		dcc.smearP3[sec][sl]    = data[row][4];
		dcc.smearP4[sec][sl]    = data[row][5];
		dcc.smearScale[sec][sl] = data[row][6];
		
		if(dcc.smearScale[sec][sl] > 1)
		{
			cout << "  !!!! DC Warning: the smearing parameter is greater than one for sector " << sec << " sl " << sl
			<< ". That means that the DC response in GEMC will have"
			<< " worse resoultion than the data. " << endl;
		}

	}

	
	// reading DC core parameters
	database   = "/geometry/dc/superlayer";
	auto_ptr<Assignment> dcCoreModel(calib->GetAssignment(database));
	for(size_t rowI = 0; rowI < dcCoreModel->GetRowsCount(); rowI++)
		dcc.dLayer[rowI] = dcCoreModel->GetValueDouble(rowI, 6);
	
	dcc.driftVelocity[0] = dcc.driftVelocity[1] = 0.053;  ///< drift velocity is 53 um/ns for region1
	dcc.driftVelocity[2] = dcc.driftVelocity[3] = 0.026;  ///< drift velocity is 26 um/ns for region2
	dcc.driftVelocity[4] = dcc.driftVelocity[5] = 0.036;  ///< drift velocity is 36 um/ns for region3
	
	for(int l=0; l<6; l++)
		dcc.miniStagger[l] = 0;


	// loading translation table
	dcc.TT = TranslationTable("dcTT");
	cout << "  > Data loaded in translation table " << dcc.TT.getName() << endl;

	// setting voltage signal parameters
	dcc.vpar[0] = 50;  // delay, ns
	dcc.vpar[1] = 10;  // rise time, ns
	dcc.vpar[2] = 20;  // fall time, ns
	dcc.vpar[3] = 1;   // amplifier


	return dcc;
}


map<string, double> dc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	int SECI  = identity[0].id - 1;
	int SLI   = identity[1].id - 1;
	int nwire = identity[3].id;
	
	// nwire position information
	double ylength =  aHit->GetDetector().dimensions[3];    ///< G4Trap Semilength
	double deltay  = 2.0*ylength/dcc.NWIRES;                ///< Y length of cell
	double WIRE_Y  = nwire*deltay + dcc.miniStagger[SLI];   ///< Center of wire hit
	
	vector<int>           stepTrackId = aHit->GetTIds();
	vector<double>        stepTime    = aHit->GetTime();
	vector<G4double>      Edep        = aHit->GetEdep();
	vector<G4ThreeVector> pos         = aHit->GetPos();
	vector<G4ThreeVector> Lpos        = aHit->GetLPos();
	
	unsigned nsteps = Edep.size();

	// Identifying the fastest - given by time + doca(s) / drift velocity
	// trackId Strong and Weak in case the track does or does not deposit enough energy
	int trackIds = -1;
 	int trackIdw = -1;
	double minTime  = 10000;
	double signal_t = 0;

	for(unsigned int s=0; s<nsteps; s++)
	{
		G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z()); // local cylinder
		signal_t = stepTime[s]/ns + DOCA.mag()/dcc.driftVelocity[SLI];

		// threshold hardcoded, please get from parameters
		if(signal_t < minTime)
		{
			trackIdw = stepTrackId[s];
			minTime = signal_t;
			
			if(Edep[s] >= dcc.dcThreshold*eV)
			{
				trackIds = stepTrackId[s];
			}
		}
	}
	
	// If no step pass the threshold, getting the fastest signal of the weak tracks
	if(trackIds == -1)
		trackIds = trackIdw;
	
	// Left / Right ambiguity
	// Finding DOCA
	double LR      = 0;
	double doca    = 10000;
	
	for(unsigned int s=0; s<nsteps; s++)
	{
		G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z());
		
		if(DOCA.mag() <= doca && stepTrackId[s] == trackIds )
		{
			doca = DOCA.mag();
			if(DOCA.y() >=0 ) LR = 1;
			else  LR = -1;
		}
	}
	
	// percentage distance from the wire
	double X = (doca/cm) / (2*dcc.dLayer[SLI]);

	// smeading doca by DOCA dependent function. This is the sigma of doca.
	double smearF = dcc.smearP1[SECI][SLI] + dcc.smearP2[SECI][SLI]/pow(dcc.smearP3[SECI][SLI] + X, 2) + dcc.smearP4[SECI][SLI]*pow(X, 8);
	
	// smearing doca with a gaussian around doca and sigma defined above
	double sdoca = fabs(CLHEP::RandGauss::shoot(doca, smearF*dcc.smearScale[SECI][SLI]));
	
	// distance-dependent efficiency as a function of doca
	double ddEff = dcc.iScale[SLI]*(dcc.P1[SLI]/pow(X*X + dcc.P2[SLI], 2) + dcc.P3[SLI]/pow( (1-X) + dcc.P4[SLI], 2));
	double random = G4UniformRand();

	int ineff = 1;
	if(random < ddEff || X > 1) ineff = -1;
		
	// recording smeared and un-smeared quantities
	dgtz["hitn"]       = hitn;
	dgtz["sector"]     = identity[0].id;
	dgtz["superlayer"] = SLI+1;
	dgtz["layer"]      = identity[2].id;
	dgtz["wire"]       = nwire;
	dgtz["tdc"]        = (int) sdoca/dcc.driftVelocity[SLI]*100;
	dgtz["LR"]         = LR;
	dgtz["doca"]       = doca;
	dgtz["sdoca"]      = sdoca;
	dgtz["time"]       = ineff *doca/dcc.driftVelocity[SLI];
	dgtz["stime"]      = ineff*sdoca/dcc.driftVelocity[SLI];
	
	return dgtz;
}

// routine to determine the wire number based on the hit position
vector<identifier>  dc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	int SLI = yid[1].id - 1;
	
	G4StepPoint   *prestep   = aStep->GetPreStepPoint();
	G4StepPoint   *poststep  = aStep->GetPostStepPoint();
	
	G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction
	G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);
	
	
	double ylength = Detector.dimensions[3];  ///< G4Trap Semilength
	double deltay  = 2.0*ylength/dcc.NWIRES;
	double loc_y   = Lxyz.y() + ylength - dcc.miniStagger[SLI];    ///< Distance from bottom of G4Trap - modified by ministaggger
	
	int nwire = (int) floor(loc_y/deltay);
	
	// resetting nwire for extreme cases
	if(nwire <= 0 )  nwire = 1;
	if(nwire == 113) nwire = 112;
	
	// setting wire number
	yid[3].id = nwire;
	
	// checking that the next wire is not the one closer
	if(fabs( (nwire+1)*deltay - loc_y ) < fabs( nwire*deltay - loc_y ) && nwire != 112 )
		yid[3].id = nwire + 1;
	
	// all energy to this wire (no energy sharing)
	yid[3].id_sharing = 1;
	
	return yid;
}

void dc_HitProcess::initWithRunNumber(int runno)
{
	if(dcc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		dcc = initializeDCConstants(runno);
		dcc.runNo = runno;
	}
}


map< string, vector <int> >  dc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > dc_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double dc_HitProcess :: voltage(double charge, double time, double forTime)
{
	//	return 0.0;
	return DGauss(forTime, dcc.vpar, charge, time);
}



// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> dc_HitProcess :: electronicNoise()
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

// this static function will be loaded first thing by the executable
dcConstants dc_HitProcess::dcc = initializeDCConstants(-1);







