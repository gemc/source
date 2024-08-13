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


static dcConstants initializeDCConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	dcConstants dcc;
	
	// used in process ID so defined here
	dcc.NWIRES = 113;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return dcc;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	dcc.runNo = runno;
	dcc.dcThreshold  = 50;  // eV
	
	// database
	if(getenv ("CCDB_CONNECTION") != nullptr)
		dcc.connection = (string) getenv("CCDB_CONNECTION");
	else
		dcc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(dcc.connection));
	
	
	// reading efficiency parameters
	snprintf(dcc.database, sizeof(dcc.database), "/calibration/dc/signal_generation/inefficiency:%d:%s%s", dcc.runNo, digiVariation.c_str(), timestamp.c_str());
	vector<vector<double> > data;
	calib->GetCalib(data, dcc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		int sec = data[row][0] - 1;
		int sl  = data[row][1] - 1;
		dcc.iScale[sec][sl] = data[row][3];
		dcc.P1[sec][sl]     = data[row][4];
		dcc.P2[sec][sl]     = data[row][5];
		dcc.P3[sec][sl]     = data[row][6];
		dcc.P4[sec][sl]     = data[row][7];
	}
	
	// reading smearing parameters
	snprintf(dcc.database, sizeof(dcc.database),  "/calibration/dc/signal_generation/doca_smearing:%d:%s%s", dcc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear();
	calib->GetCalib(data, dcc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		int sec = data[row][0] - 1;
		int sl  = data[row][1] - 1;
		dcc.smearP0[sec][sl]    = data[row][3];
		dcc.smearP1[sec][sl]    = data[row][4];
		dcc.smearP2[sec][sl]    = data[row][5];
		dcc.smearP3[sec][sl]    = data[row][6];
		dcc.smearP4[sec][sl]    = data[row][7];
	}
	
	
	//********************************************
	//calculating distance to time:
	snprintf(dcc.database, sizeof(dcc.database),  "/calibration/dc/time_to_distance/time2dist:%d:%s%s", dcc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear();
	calib->GetCalib(data, dcc.database);
	
	for(unsigned row = 0; row < data.size(); row++) {
		int sec = data[row][0] - 1;
		int sl  = data[row][1] - 1;
		dcc.v0[sec][sl] = data[row][3];
		dcc.deltanm[sec][sl] = data[row][4]; //used in exponential function only
		dcc.tmaxsuperlayer[sec][sl] = data[row][5];
		// Row left out, corresponds to distbfield
		dcc.delta_bfield_coefficient[sec][sl] = data[row][7];
		dcc.deltatime_bfield_par1[sec][sl] = data[row][8];
		dcc.deltatime_bfield_par2[sec][sl] = data[row][9];
		dcc.deltatime_bfield_par3[sec][sl] = data[row][10];
		dcc.deltatime_bfield_par4[sec][sl] = data[row][11];
		dcc.R[sec][sl] = data[row][13];     // used in polynomial function only
		dcc.vmid[sec][sl] = data[row][14];  // used in polynomial function only
		//        cout << dcc.v0[sec][sl] << " " << dcc.deltanm[sec][sl] << " " << dcc.tmaxsuperlayer[sec][sl] << " " << dcc.R[sec][sl] << " " << dcc.vmid[sec][sl] << endl;
	}
	
	
	
	// T0 corrections: a delay to be introduced (plus sign) to the TDC timing
	snprintf(dcc.database, sizeof(dcc.database),  "/calibration/dc/time_corrections/T0Corrections:%d:%s%s", dcc.runNo, digiVariation.c_str(), timestamp.c_str());
	data.clear();
	calib->GetCalib(data,  dcc.database);
	
	for(unsigned row = 0; row < data.size(); row++) {
		int sec    = data[row][0] - 1;
		int sl     = data[row][1] - 1;
		int slot   = data[row][2] - 1;
		int cable  = data[row][3] - 1;
		dcc.T0Correction[sec][sl][slot][cable] = data[row][4];
	}
	//********************************************
	
	
	
	// reading DC core parameters
	snprintf(dcc.database, sizeof(dcc.database),  "/geometry/dc/superlayer:%d:%s%s", dcc.runNo, digiVariation.c_str(), timestamp.c_str());
	unique_ptr<Assignment> dcCoreModel(calib->GetAssignment(dcc.database));
	for(size_t rowI = 0; rowI < dcCoreModel->GetRowsCount(); rowI++){
		dcc.dLayer[rowI] = dcCoreModel->GetValueDouble(rowI, 6);
	}
	
	dcc.dmaxsuperlayer[0] = 2*dcc.dLayer[0];
	dcc.dmaxsuperlayer[1] = 2*dcc.dLayer[1];
	dcc.dmaxsuperlayer[2] = 2*dcc.dLayer[2];
	dcc.dmaxsuperlayer[3] = 2*dcc.dLayer[3];
	dcc.dmaxsuperlayer[4] = 2*dcc.dLayer[4];
	dcc.dmaxsuperlayer[5] = 2*dcc.dLayer[5];
	
	
	// even number layers are closer to the beamline:
	// layers 1,3,5 have +300 micron
	// layers 2,4,6 have -300 micron
	dcc.miniStagger[0] =  0.300*mm;
	dcc.miniStagger[1] = -0.300*mm;
	dcc.miniStagger[2] =  0.300*mm;
	dcc.miniStagger[3] = -0.300*mm;
	dcc.miniStagger[4] =  0.300*mm;
	dcc.miniStagger[5] = -0.300*mm;
	
	
	// loading translation table; CURRENTLY NOT USED
	dcc.TT = TranslationTable("dcTT");
	cout << "  > Data loaded in translation table " << dcc.TT.getName() << endl;
	
	// setting voltage signal parameters; CURRENTLY NOT USED
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
	rejectHitConditions = false;
	writeHit = true;
	
	if(aHit->isBackgroundHit == 1) {
		
		vector<double>        stepTime    = aHit->GetTime();
		//	cout << " This is a background hit with time " << stepTime[0] << endl;
		
		//		int SECI  = identity[0].id - 1;
		int SLI   = identity[1].id - 1;
		//		int LAY   = identity[2].id - 1;
		int nwire = identity[3].id;
		
		dgtz["hitn"]       = hitn;
		dgtz["sector"]     = identity[0].id;
		dgtz["layer"]      = SLI*6 + identity[2].id;
		dgtz["component"]  = nwire;
		dgtz["TDC_order"] = 0;
		dgtz["TDC_TDC"]   = stepTime[0];
		
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
	double hit_signal_t = 0;
	
	for(unsigned int s=0; s<nsteps; s++)
	{
		G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z()); // local cylinder
		signal_t = stepTime[s]/ns + DOCA.mag()/(dcc.v0[SECI][SLI]*cm/ns);
		// cout << "signal_t: " << signal_t << " stepTime: " << stepTime[s] << " DOCA: " << DOCA.mag() << " driftVelocity: " << dcc.driftVelocity[SLI] << " Lposy: " << Lpos[s].y() << " ylength: " << ylength << " WIRE_Y: " << WIRE_Y << " Lposz: " << Lpos[s].z() << " dcc.NWIRES: " << dcc.NWIRES << endl;
		
		if(signal_t < minTime)
		{
			trackIdw = stepTrackId[s];
			minTime = signal_t;
			
			if(Edep[s] >= dcc.dcThreshold*eV) {
				// new hit time
				// (w/o the drift time)
				// not activated yet
				// hit_signal_t = stepTime[s]/ns;
				
				trackIds = stepTrackId[s];
			}
		}
	}
	
	
	
	// If no step pass the threshold, getting the fastest signal with no threshold: FOR MAC, IS THIS WHAT WE WANT?
	if(trackIds == -1)
		trackIds = trackIdw;
	
	
	// Left / Right ambiguity (notice: not used anymore?)
	// Finding DOCA
	// double LR       = 0;
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
			
			//B-field correction: correct alpha with theta0, the angle corresponding to the isochrone lines twist due to the electric field
			thisMgnf = mgnf[s]/tesla; // Given in Tesla
			double theta0 = acos(1-0.02*thisMgnf)/deg;
			alpha-= dcc.fieldPolarity*theta0;
			
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
			// LR not used anymore
			// if(DOCA.y() >=0 ) LR = 1;
			// else  LR = -1;
			
			//Get beta-value of the particle:
			beta_particle = mom[s].mag()/E[s];
			
		}
	}
	//******************************************************************
	
	// percentage distance from the wire
	double X = (doca/cm) / (2*dcc.dLayer[SLI]);
	
	// distance-dependent fractional inefficiency as a function of doca
	double ddEff = dcc.iScale[SECI][SLI]*(dcc.P1[SECI][SLI]/pow(X*X + dcc.P2[SECI][SLI], 2) + dcc.P3[SECI][SLI]/pow( (1-X) + dcc.P4[SECI][SLI], 2));
	double random = G4UniformRand();
	
	// unsmeared time, based on the dist-time-function and alpha;
	double unsmeared_time = calc_Time(doca/cm,dcc.dmaxsuperlayer[SLI],dcc.tmaxsuperlayer[SECI][SLI],alpha,thisMgnf,SECI,SLI);
	
	// Include time smearing calculated from doca resolution
	double dt_random_in = doca_smearing(X, beta_particle, SECI, SLI);
	//double dt_random = dt_random_in*CLHEP::RandLandau::shoot();
	double dt_random = fabs(CLHEP::RandGauss::shoot(0,dt_random_in));
	//cout << X << " " << beta_particle << " " << dcc.v0[SECI][SLI] << " " << dt_random_in << endl;
	
	// Now calculate the smeared time:
	// adding the time of hit from the start of the event (signal_t), which also has the drift velocity into it
	double smeared_time = unsmeared_time + dt_random + hit_signal_t + dcc.get_T0(SECI, SLI, LAY, nwire);
	
	// cout << " DC TIME stime: " << smeared_time << " X: " << X << "  doca: " << doca/cm << "  dmax: " << dcc.dmaxsuperlayer[SLI] << "    tmax: " << dcc.tmaxsuperlayer[SECI][SLI] << "   alpha: " << alpha << "   thisMgnf: " << thisMgnf << " SECI: " << SECI << " SLI: " << SLI << endl;
	
	int ineff = 1;
	if(random < ddEff || X > 1) ineff = -1;
	
	// recording smeared and un-smeared quantities
	dgtz["hitn"]       = hitn;
	dgtz["sector"]     = identity[0].id;
	dgtz["layer"]      = SLI*6 + identity[2].id;
	dgtz["component"]  = nwire;
	dgtz["TDC_order"]  = 0;
	dgtz["TDC_TDC"]    = smeared_time;
	
	// decide if write an hit or not based on inefficiency value
	rejectHitConditions=(ineff==-1);

	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}

	// cout << " ASD hit n " << hitn << " reject: " << rejectHitConditions << " ddeff:" << ddEff << " random:" << random << " X:" << X << " writeHit: " << writeHit << endl;

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



map< string, vector <int> >  dc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

// OLD Exponential function: returns a time in ns give:
// x      = distance from the wire, in cm
// dmax   = cell size in superlayer
// tmax   = t max in superlayer
// alpha  = local angle of the track
// bfield = magnitude of field in tesla
// sector      = sector
// superlayer      = superlayer
double dc_HitProcess :: calc_Time_exp(double x, double dmax, double tmax, double alpha, double bfield, int sector, int superlayer)
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

// NEW Polynomial function: returns a time in ns give:
// x      = distance from the wire, in cm
// dmax   = cell size in superlayer
// tmax   = t max in superlayer
// alpha  = local angle of the track
// bfield = magnitude of field in tesla
// sector      = sector
// superlayer      = superlayer
double dc_HitProcess :: calc_Time(double x, double dmax, double tmax, double alpha, double bfield, int sector, int superlayer)
{
	if(x>dmax)
		x=dmax;
	double time = 0;
	// alpha correction
	double cos30minusalpha=(double)cos((30. - alpha)*deg);
	double dmaxalpha = dmax*cos30minusalpha;
	double xhatalpha = x/dmaxalpha;
	//   rcapital is an intermediate parameter
	double rcapital = dcc.R[sector][superlayer]*dmax;
	//   delt is another intermediate parameter
	double delt=tmax-dmax/dcc.v0[sector][superlayer];
	double delv=1./dcc.vmid[sector][superlayer]-1./dcc.v0[sector][superlayer];
	//   now calculate the primary parameters a, b, c, d
	
	double c = ((3.*delv)/(dcc.R[sector][superlayer]*dmax)+
					(12*dcc.R[sector][superlayer]*dcc.R[sector][superlayer]*delt)/
					(2.*(1-2*dcc.R[sector][superlayer])*(dmax*dmax)));
	c = c /(4.-(1.-6.*dcc.R[sector][superlayer]*dcc.R[sector][superlayer])/(1.-2.*dcc.R[sector][superlayer]));
	double b = delv/(rcapital*rcapital) - 4.*c/(3.*rcapital);
	double d = 1./dcc.v0[sector][superlayer];
	double a = (tmax -  b*dmaxalpha*dmaxalpha*dmaxalpha -
					c*dmaxalpha*dmaxalpha - d*dmaxalpha)/(dmaxalpha*dmaxalpha*dmaxalpha*dmaxalpha) ;
	time = a*x*x*x*x + b*x*x*x + c*x*x + d*x ;
	
	double deltatime_bfield = dcc.delta_bfield_coefficient[sector][superlayer]*pow(bfield,2)*tmax*(dcc.deltatime_bfield_par1[sector][superlayer]*xhatalpha+dcc.deltatime_bfield_par2[sector][superlayer]*pow(xhatalpha, 2)+ dcc.deltatime_bfield_par3[sector][superlayer]*pow(xhatalpha, 3)+dcc.deltatime_bfield_par4[sector][superlayer]*pow(xhatalpha, 4));
	
	//cout<<"dt = "<<deltatime_bfield<<" C0 "<<dcc.delta_bfield_coefficient[sector][superlayer]<<endl;
	//calculate the time at alpha deg. and at a non-zero bfield
	time += deltatime_bfield;
	return time;
}

// Define DOCA smearing based on data parameterization
// x: distance from the wire normalized to the cell size
// beta: beta of the particle
// sector: DC sector
// superlayer: DC superlayer
// returns time smearing in ns
double dc_HitProcess :: doca_smearing(double x, double beta, int sector, int superlayer){
	double doca_smear = 0.0;
	
	double dmax = 1;
	if(x>dmax) x=dmax;
	
	doca_smear  = (dcc.smearP0[sector][superlayer]
						+ dcc.smearP1[sector][superlayer] * x
						+ dcc.smearP2[sector][superlayer] * x * x
						+ dcc.smearP3[sector][superlayer] * x * x * x
						+ dcc.smearP4[sector][superlayer] * x * x * x * x) * cm;
	doca_smear  = doca_smear/(dcc.v0[sector][superlayer]*cm/ns);
	
	return doca_smear;
}


// - electronicNoise: returns a vector of hits generated / by electronics: NOT CURRENTLY USED
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


void dc_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	string hardcodedAsciiTorusMapName  = "TorusSymmetric";
	string hardcodedBinaryTorusMapName = "binary_torus";

	if(dcc.runNo != runno) {
		
		dcc = initializeDCConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		dcc.runNo = runno;
		
		double scaleFactor = 1;
		vector<aopt> FIELD_SCALES_OPTION = gemcOpt.getArgs("SCALE_FIELD");
		for (unsigned int f = 0; f < FIELD_SCALES_OPTION.size(); f++) {
			vector < string > scales = getStringVectorFromStringWithDelimiter(FIELD_SCALES_OPTION[f].args, ",");
			if(scales.size() == 2) {
				if (scales[0].find(hardcodedAsciiTorusMapName) != string::npos) {
					scaleFactor = get_number(scales[1]);
				} else 	if (scales[0].find(hardcodedBinaryTorusMapName) != string::npos) {
					scaleFactor = get_number(scales[1]);
				}
			}
		}
		dcc.fieldPolarity = scaleFactor;
		
		string nofield = gemcOpt.optMap["NO_FIELD"].args;
		if(nofield == "all") {
			dcc.fieldPolarity = 1;
		}
		
		cout << " > Initializing " << HCname << " digitization for run number: " << dcc.runNo << ", torus polarity: " << dcc.fieldPolarity << endl;
		
		
	}
}

// this static function will be loaded first thing by the executable
dcConstants dc_HitProcess::dcc = initializeDCConstants(-1);
