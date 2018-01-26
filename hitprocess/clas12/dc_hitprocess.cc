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
    //	dcc.date       = "2016-03-15";
    dcc.date       = "2017-08-01";
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
    
    
    //********************************************
    //calculating distance to time:
    database  = "/calibration/dc/time_to_distance/tvsx_devel_v2";
    data.clear();
    calib->GetCalib(data, database);
    
    for(unsigned row = 0; row < data.size(); row++)
    {
        int sec = data[row][0] - 1;
        int sl  = data[row][1] - 1;
        dcc.v0[sec][sl] = data[row][2];
        dcc.deltanm[sec][sl] = data[row][3];
        dcc.tmaxsuperlayer[sec][sl] = data[row][4];
        //Row 5 is left out, corresponds to distbfield
        dcc.delta_bfield_coefficient[sec][sl] = data[row][6];
        dcc.deltatime_bfield_par1[sec][sl] = data[row][7];
        dcc.deltatime_bfield_par2[sec][sl] = data[row][8];
        dcc.deltatime_bfield_par3[sec][sl] = data[row][9];
        dcc.deltatime_bfield_par4[sec][sl] = data[row][10];
    }
    
    
    //Include smearing parameters for time walks: (up to now fixed values, will be included in ccdb soon):
    dcc.smear_time_walk[0] = 0.001; //Corresponds to smearing at low distances from the wire (given in mm)
    dcc.smear_time_walk[1] = 3.0; //Adjusts ratio between time walk effects close and far away from the wire
    dcc.smear_time_walk[2] = 30; //Adjusts random time walk smearing
    //********************************************
    
    
    
    // reading DC core parameters
    database   = "/geometry/dc/superlayer";
    auto_ptr<Assignment> dcCoreModel(calib->GetAssignment(database));
    for(size_t rowI = 0; rowI < dcCoreModel->GetRowsCount(); rowI++){
        dcc.dLayer[rowI] = dcCoreModel->GetValueDouble(rowI, 6);
        dcc.driftVelocity[rowI] = dcCoreModel->GetValueDouble(rowI, 7);
    }
   
    dcc.dmaxsuperlayer[0] = 2*dcc.dLayer[0];
    dcc.dmaxsuperlayer[1] = 2*dcc.dLayer[1];
    dcc.dmaxsuperlayer[2] = 2*dcc.dLayer[2];
    dcc.dmaxsuperlayer[3] = 2*dcc.dLayer[3];
    dcc.dmaxsuperlayer[4] = 2*dcc.dLayer[4];
    dcc.dmaxsuperlayer[5] = 2*dcc.dLayer[5];

    dcc.driftVelocity[0] = dcc.driftVelocity[1] = 0.053;  ///< drift velocity is 53 um/ns for region1
	dcc.driftVelocity[2] = dcc.driftVelocity[3] = 0.026;  ///< drift velocity is 26 um/ns for region2
	dcc.driftVelocity[4] = dcc.driftVelocity[5] = 0.036;  ///< drift velocity is 36 um/ns for region3

    
    // even closer:
    // layers 1,3,5 have +300 micron
    // layers 2,4,6 have -300 micron
    dcc.miniStagger[0] =  0.300*mm;
    dcc.miniStagger[1] = -0.300*mm;
    dcc.miniStagger[2] =  0.300*mm;
    dcc.miniStagger[3] = -0.300*mm;
    dcc.miniStagger[4] =  0.300*mm;
    dcc.miniStagger[5] = -0.300*mm;
    
    
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

	if(aHit->isBackgroundHit == 1) {

		vector<double>        stepTime    = aHit->GetTime();
	//	cout << " This is a background hit with time " << stepTime[0] << endl;
	// recording smeared and un-smeared quantities
		dgtz["hitn"]       = hitn;
		dgtz["sector"]     = identity[0].id;
		dgtz["layer"]      = identity[1].id;
		dgtz["wire"]       = identity[2].id;;
		dgtz["tdc"]        = stepTime[0];
		return dgtz;
	}


    int SECI  = identity[0].id - 1;
    int SLI   = identity[1].id - 1;
    int LAY   = identity[2].id - 1;
    int nwire = identity[3].id;
    
    // nwire position information
    double ylength =  aHit->GetDetector().dimensions[3];    ///< G4Trap Semilength
    double deltay  = 2.0*ylength/dcc.NWIRES;                ///< Y length of cell
    double WIRE_Y  = nwire*deltay;                          ///< Center of wire hit
    if(SLI > 3) WIRE_Y += dcc.miniStagger[LAY];             ///< Region 3 (SLI 4 and 5) have mini-stagger for the sense wires
    
    vector<int>           stepTrackId = aHit->GetTIds();
    vector<double>        stepTime    = aHit->GetTime();
    vector<double>        mgnf        = aHit->GetMgnf();
    vector<G4double>      Edep        = aHit->GetEdep();
    vector<G4ThreeVector> pos         = aHit->GetPos();
    vector<G4ThreeVector> Lpos        = aHit->GetLPos();
    vector<G4ThreeVector> mom         = aHit->GetMoms();
    vector<double>        E           = aHit->GetEs();
    
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
            
            if(Edep[s] >= dcc.dcThreshold*eV) {
                trackIds = stepTrackId[s];
            }
        }
    }
    
    
    
    // If no step pass the threshold, getting the fastest signal of the weak tracks
    if(trackIds == -1)
        trackIds = trackIdw;
    
    // Left / Right ambiguity
    // Finding DOCA
    double LR       = 0;
    double doca     = 10000;
    double thisMgnf = 0;
    
    //Quantities needed for the calculation of alpha:
    //******************************************************************
    double alpha = 0;
    double rotate_to_sec = (-60*SECI)*deg; //Rotation towards fireing sector
    
    int sl_sign; //check superlayer and define sign --> Needed for orientation within the wire system
    for(int i=0;i<3;i++){
        if(SLI == 2*i+1){
            sl_sign = -1;
        }else sl_sign = 1;
    }
    
    double rotate_to_wire = 6*sl_sign*deg; //Angle for the final rotation into the wire system
    double const2 = sl_sign*sin(6*deg);
    double const1= cos(6*deg);
    double beta_particle = 0.0; //beta-value of the particle --> Needed for determining time walks
    
    G4ThreeVector rotated_vector;
    
    for(unsigned int s=0; s<nsteps; s++)
    {
        G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z());
        
        if(DOCA.mag() <= doca && stepTrackId[s] == trackIds )
        {
            //Get the momentum vector, which shall be rotated:
            rotated_vector = G4ThreeVector(mom[s].x(),mom[s].y(),mom[s].z());
            
            //First, rotate px and py to the sector,that has fired:
            rotated_vector.rotateZ(rotate_to_sec);
            
            //Secondly, rotate px' and pz towards the firing sector:
            rotated_vector.rotateY(-25*deg);
            
            //Finally, rotate px'' and py' into the wire coordinate system:
            rotated_vector.rotateZ(rotate_to_wire);
            
            //Now calculate alpha according to Macs definition:
            alpha = asin((const1*rotated_vector.x() + const2*rotated_vector.y())/rotated_vector.mag())/deg;
            
            
            // compute reduced alpha (VZ)
            // alpha in radians
            double ralpha = fabs(alpha*deg);
            
            while (ralpha > pi / 3.) {
                ralpha -= pi / 3.;
            }
            if (ralpha > pi / 6.) {
                ralpha = pi / 3. - ralpha;
            }
            //alpha in degrees (reduced alpha always between 0 and 30 deg.)
            alpha = ralpha/deg;
            
            doca = DOCA.mag();
            if(DOCA.y() >=0 ) LR = 1;
            else  LR = -1;
            thisMgnf = mgnf[s]; //Given in Tesla
            
            //Get beta-value of the particle:
            beta_particle = mom[s].mag()/E[s];
            
        }
    }
    //******************************************************************
    
    // percentage distance from the wire
    double X = (doca/cm) / (2*dcc.dLayer[SLI]);
    
    // distance-dependent efficiency as a function of doca
    double ddEff = dcc.iScale[SLI]*(dcc.P1[SLI]/pow(X*X + dcc.P2[SLI], 2) + dcc.P3[SLI]/pow( (1-X) + dcc.P4[SLI], 2));
    double random = G4UniformRand();
    
    //unsmeared time, based on the dist-time-function and alpha;
    double unsmeared_time = calc_Time(doca/cm,dcc.dmaxsuperlayer[SLI],dcc.tmaxsuperlayer[SECI][SLI],alpha,thisMgnf/tesla,SECI,SLI);
    
    // Now include (random) time walk contributions:
    
    // Include ionisation effects:
    double dt_walk_in = time_walk_core(doca/cm,dcc.dmaxsuperlayer[SLI],dcc.smear_time_walk[0]*beta_particle*beta_particle,dcc.smear_time_walk[1],1,dcc.v0[SECI][SLI]);
    //Translate average value into detector response/reported value (landau-function):
    double dt_walk = dt_walk_in + 0.5*dt_walk_in*CLHEP::RandLandau::shoot();
    
    // Include intrinsic random time walks (due to multiple scattering):
    double dt_random_in = time_rnd_core(doca/cm,dcc.smear_time_walk[2],dcc.v0[SECI][SLI]);
    //And translate it to a proper detector response/reported value (gaussian):
    double dt_random = CLHEP::RandGauss::shoot(0,dt_random_in);
    
    //Now calculate the smeared time:
    double smeared_time = unsmeared_time + dt_walk + dt_random;
    
    int ineff = 1;
    if(random < ddEff || X > 1) ineff = -1;
    
    // recording smeared and un-smeared quantities
    dgtz["hitn"]       = hitn;
    dgtz["sector"]     = identity[0].id;
    dgtz["layer"]      = SLI*6 + identity[2].id;
    dgtz["wire"]       = nwire;
    dgtz["tdc"]        = smeared_time;
    dgtz["LR"]         = LR;
    dgtz["doca"]       = doca;
    dgtz["time"]       = ineff*unsmeared_time;
    dgtz["stime"]      = ineff*smeared_time;
    
    // decide if write an hit or not
    writeHit = true;
    // define conditions to reject hit
    bool rejectHitConditions = false;
    if(rejectHitConditions) {
         writeHit = false;
    }

    return dgtz;
}

// routine to determine the wire number based on the hit position
vector<identifier>  dc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
    vector<identifier> yid = id;
    
    G4StepPoint   *prestep   = aStep->GetPreStepPoint();
    G4StepPoint   *poststep  = aStep->GetPostStepPoint();
    G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction
    G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
    ->GetTopTransform().TransformPoint(xyz);
    
    
    
    
    double ylength = Detector.dimensions[3];  ///< G4Trap Semilength
    double deltay  = 2.0*ylength/dcc.NWIRES;
    double loc_y   = Lxyz.y() + ylength;    ///< Distance from bottom of G4Trap. ministaggger does not affect it since the field/guardwires are fixed.
    
    int nwire = (int) floor(loc_y/deltay);
    
    // resetting nwire for extreme cases
    if(nwire <= 0 )  nwire = 1;
    if(nwire >= 113) nwire = 112;
    
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

// returns a time in ns give:
// x      = distance from the wire, in cm
// dmax   = cell size in superlayer
// tmax   = t max in superlayer
// alpha  = polar angle of the track
// bfield = magnitude of field in tesla
// sector      = sector
// superlayer      = superlayer
double dc_HitProcess :: calc_Time(double x, double dmax, double tmax, double alpha, double bfield, int sector, int superlayer)
{
    
    double rtime = 0.0;
    double FracDmaxAtMinVel = 0.615;
    // Assume a functional form (time=x/v0+a*(x/dmax)**n+b*(x/dmax)**m)
    // for time as a function of x for theta = 30 deg.
    // first, calculate n
    double n = ( 1.+ (dcc.deltanm[sector][superlayer]-1.)*pow(FracDmaxAtMinVel, dcc.deltanm[sector][superlayer]) )/( 1.- pow(FracDmaxAtMinVel, dcc.deltanm[sector][superlayer]));
    //now, calculate m
    double m = n + dcc.deltanm[sector][superlayer];
    // determine b from the requirement that the time = tmax at dist=dmax
    double b = (tmax - dmax/dcc.v0[sector][superlayer])/(1.- m/n);
    
    // determine a from the requirement that the derivative at
    // d=dmax equal the derivative at d=0
    double a = -b*m/n;
    
    double cos30minusalpha=(double)cos((30. - alpha)*deg);
    double xhat = x/dmax;
    double dmaxalpha = dmax*cos30minusalpha;
    double xhatalpha = x/dmaxalpha;
    
    //     now calculate the dist to time function for theta = 'alpha' deg.
    //     Assume a functional form with the SAME POWERS N and M and
    //     coefficient a but a new coefficient 'balpha' to replace b.
    //     Calculate balpha from the constraint that the value
    //     of the function at dmax*cos30minusalpha is equal to tmax
    
    //     parameter balpha (function of the 30 degree paramters a,n,m)
    double balpha = ( tmax - dmaxalpha/dcc.v0[sector][superlayer] - a*pow(cos30minusalpha,n))/pow(cos30minusalpha, m);
    
    //      now calculate function
    rtime = x/dcc.v0[sector][superlayer] + a*pow(xhat, n) + balpha*pow(xhat, m);
    
    double deltatime_bfield = dcc.delta_bfield_coefficient[sector][superlayer]*pow(bfield,2)*tmax*(dcc.deltatime_bfield_par1[sector][superlayer]*xhatalpha+dcc.deltatime_bfield_par2[sector][superlayer]*pow(xhatalpha, 2)+ dcc.deltatime_bfield_par3[sector][superlayer]*pow(xhatalpha, 3)+dcc.deltatime_bfield_par4[sector][superlayer]*pow(xhatalpha, 4));
    
    //cout<<"dt = "<<deltatime_bfield<<" C0 "<<dcc.delta_bfield_coefficient[sector][superlayer]<<endl;
    //calculate the time at alpha deg. and at a non-zero bfield
    rtime += deltatime_bfield;
    
    return rtime;
}

//Taking care of time walks due to discrete ionization processes:
// x       = distance from the wire, in cm
// dmax    = cell size in superlayer
// epsilon = adjustable factor (in mm) times beta² of the particle, the factor is adjusted according to data at small distances from the wire
// R       = Another (relative) parameter to adjust data at large distances from the wire
// kappa   = Parameter to adjust the matching between the two time walk distributions --> Should not be touched!
// v0      = velocity of the particle
//The time returned is given in ns
double dc_HitProcess :: time_walk_core(double x, double dmax, double epsilon, double R, double kappa, double v0){
    double out_walk = 0.0;
    
    if(epsilon > 0 && v0 > 0){//We dont want any trouble with 1/0...
        double xcrit = 0.615*dmax; //Reflection point
        //Two extreme cases have to be considered:
        
        //i) x < xcrit: Distances close to the wire
        double dt_close = (sqrt(epsilon*epsilon + x*x) - x)/v0; //This expression is based on geometrical considerations (i.e. ions along a particle track)
        
        //ii) x >= xcrit: Distances far from the wire
        double dt_far = R*epsilon*epsilon/(v0*((dmax-x) + epsilon)); //Using the approach electric field ~ 1/r
        
        //Now merge both expression via a sigmoid function, in order to have a continious distribution without discrete steps:
        double arg = 1 + exp(-kappa*(x-xcrit));
        out_walk = (dt_far - dt_close)/arg + dt_close;
    }
    return out_walk;
}

//Include random walk contributions (basically scattering):
// x      = distance from the wire, in cm
// f      = factor (in 10⁻3) to adjust the random contributions to data
//v0      = velocity of the particle
//The time returned is given in ns
double dc_HitProcess :: time_rnd_core(double x, double f, double v0){
    double out_rnd = 0.0;
    
    if(x >= 0 && v0 > 0){
        out_rnd = f*(1e-3)*sqrt(x)/v0;
    }
    return out_rnd;
}

//Unused functions/objects --> Will be removed with gemc 3.0
//****************************************************************************************

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
dcConstants dc_HitProcess::dcc = initializeDCConstants(1);

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
G4ThreeVector dc_HitProcess :: psmear(G4ThreeVector p)
{
    G4ThreeVector y(p);
    y.setX(p.x() + 1);
    
    return y;
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
   return DGauss(forTime, dcc.vpar, charge, time);
}

//****************************************************************************************






