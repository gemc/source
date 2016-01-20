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
	
	// database
	dcc.runNo = runno;
	dcc.date       = "2015-11-15";
	dcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	dcc.variation  = "main";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(dcc.connection));

	
	// reading efficiency function
	string database   = "/calibration/drift_chamber/distance_dependent_inefficiency";
	auto_ptr<Assignment> ddiModel(calib->GetAssignment(database));
	// third column of this table is the eff pars
	dcc.P1 = ddiModel->GetValueDouble(0, 2);
	dcc.P2 = ddiModel->GetValueDouble(1, 2);
	dcc.P3 = ddiModel->GetValueDouble(2, 2);
	dcc.P4 = ddiModel->GetValueDouble(3, 2);
	
	// reading DC core parameters
	database   = "/geometry/dc/superlayer";
	auto_ptr<Assignment> dcCoreModel(calib->GetAssignment(database));
	for(size_t rowI = 0; rowI < dcCoreModel->GetRowsCount(); rowI++)
		dcc.dLayer[rowI] = dcCoreModel->GetValueDouble(rowI, 6);

	
	dcc.NWIRES = 113;
	
	dcc.docaSmearing = 0.3; // mm
	dcc.dcThreshold  = 50;  // eV
	dcc.driftVelocity[0] = dcc.driftVelocity[1] = 0.053;  ///< drift velocity is 53 um/ns for region1
	dcc.driftVelocity[2] = dcc.driftVelocity[3] = 0.026;  ///< drift velocity is 26 um/ns for region2
	dcc.driftVelocity[4] = dcc.driftVelocity[5] = 0.036;  ///< drift velocity is 36 um/ns for region3
	
	for(int l=0; l<6; l++)
		dcc.miniStagger[l] = 0;
	
	return dcc;
}


map<string, double> dc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	int sector   = identity[0].id;
	int SL       = identity[1].id;
	int LayerI   = identity[2].id - 1; // layer index
	int nwire    = identity[3].id;
	
	
	// nwire position information
	double ylength =  aHit->GetDetector().dimensions[3];       ///< G4Trap Semilength
	double deltay  = 2.0*ylength/dcc.NWIRES;                   ///< Y length of cell
	double WIRE_Y  = nwire*deltay + dcc.miniStagger[LayerI];   ///< Center of wire hit
	
	vector<int>           stepTrackId = aHit->GetTIds();
	vector<double>        stepTime    = aHit->GetTime();
	vector<G4double>      Edep        = aHit->GetEdep();
	vector<G4ThreeVector> pos         = aHit->GetPos();
	vector<G4ThreeVector> Lpos        = aHit->GetLPos();
	
	unsigned nsteps = Edep.size();

	// Identifying the fastest - given by time + doca(s) / drift velocity
	// trackId Strong and Weak in case the track does or does not deposit enough energy
	int trackIds, trackIdw  = -1;
	double minTime  = 10000;
	double signal_t = 0;

	// calculating doca
	for(unsigned int s=0; s<nsteps; s++)
	{
		G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z()); // local cylinder
		signal_t = stepTime[s]/ns + DOCA.mag()/dcc.driftVelocity[LayerI];
		
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
	
	// smeading doca
	double sdoca = fabs(CLHEP::RandGauss::shoot(doca, dcc.docaSmearing));  ///< smeared by 300 microns for now
	
	
	// distance-dependent efficiency as a function of doca
	
	double X = (doca/cm) / (2*dcc.dLayer[SL-1]);
	double ddEff = dcc.P1/pow(X*X + dcc.P2, 2) + dcc.P3/pow( (1-X) + dcc.P4, 2);
	double random = G4UniformRand();

	double fired = 1;
	if(random < ddEff || X > 1) fired = 0;
	
	// recording smeared and un-smeared quantities
	dgtz["hitn"]       = hitn;
	dgtz["sector"]     = sector;
	dgtz["superlayer"] = SL;
	dgtz["layer"]      = LayerI+1;
	dgtz["wire"]       = nwire;
	dgtz["LR"]         = LR;
	dgtz["doca"]       = doca;
	dgtz["sdoca"]      = sdoca;
	dgtz["time"]       =  doca/dcc.driftVelocity[LayerI];
	dgtz["stime"]      = sdoca/dcc.driftVelocity[LayerI];
	dgtz["fired"]      = fired;
	
	return dgtz;
}

// routine to determine the wire number based on the hit position
vector<identifier>  dc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	int Layer = yid[2].id;
	
	G4StepPoint   *prestep   = aStep->GetPreStepPoint();
	G4StepPoint   *poststep  = aStep->GetPostStepPoint();
	
	G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction
	G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);
	
	
	double ylength = Detector.dimensions[3];  ///< G4Trap Semilength
	double deltay  = 2.0*ylength/dcc.NWIRES;
	double loc_y   = Lxyz.y() + ylength - dcc.miniStagger[Layer-1];    ///< Distance from bottom of G4Trap - modified by ministaggger
	
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
	if(dcc.runNo != runno)
	{
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

// this static function will be loaded first thing by the executable
dcConstants dc_HitProcess::dcc = initializeDCConstants(-1);







