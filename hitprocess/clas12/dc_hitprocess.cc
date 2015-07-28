// gemc headers
#include "dc_hitprocess.h"
#include "CLHEP/Random/RandGaussT.h"


// NOTE:
// Need to run gemc with -USE_PHYSICSL=gemc
// otherwise the step size inside the cell is too small
// to calculate DOCA correctly


map<string, double> dc_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	double mini_stagger_shift;
	
	int sector   = identity[0].id;
	int SL       = identity[1].id;
	int Layer    = identity[2].id;
	int nwire    = identity[3].id;
	
	mini_stagger_shift = 0;
	if(SL > 2 && SL <= 4)     // shift only for region 2
	{
		mini_stagger_shift = mini_stagger_shift_R2;
		if(Layer == 2 || Layer == 4 || Layer == 6)  mini_stagger_shift *= -1;
	}
	else if(SL > 4)           // shift only for region 3
	{
		mini_stagger_shift = mini_stagger_shift_R3;
		if(Layer == 2 || Layer == 4 || Layer == 6)  mini_stagger_shift *= -1;
	}
	
	// nwire position information
	double ylength =  aHit->GetDetector().dimensions[3];  ///< G4Trap Semilength
	double deltay  = 2.0*ylength/NWIRES;                  ///< Y length of cell
	double WIRE_Y  = nwire*deltay + mini_stagger_shift;   ///< Center of wire hit
	
	// drift velocities
	double drift_velocity = 0;
	if(SL == 1 || SL == 2) drift_velocity = 0.053;  ///< drift velocity is 53 um/ns for region1
	if(SL == 3 || SL == 4) drift_velocity = 0.026;  ///< drift velocity is 26 um/ns for region2
	if(SL == 5 || SL == 6) drift_velocity = 0.036;  ///< drift velocity is 36 um/ns for region3
	
	// Identifying the fastest - given by time + doca(s) / drift velocity
	int trackId     = -1;
	double minTime  = 10000;
	double signal_t = 0;
	
	vector<int>           stepTrackId = aHit->GetTIds();
	vector<double>        stepTime    = aHit->GetTime();
	vector<G4double>      Edep        = aHit->GetEdep();
	vector<G4ThreeVector> pos         = aHit->GetPos();
	vector<G4ThreeVector> Lpos        = aHit->GetLPos();
	
	unsigned nsteps = Edep.size();

	// calculating doca
	for(unsigned int s=0; s<nsteps; s++)
	{
		G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z()); // local cylinder
		signal_t = stepTime[s]/ns + DOCA.mag()/drift_velocity;
		
		// threshold hardcoded, please get from parameters
		if(signal_t < minTime && Edep[s] >= 50*eV)
		{
			trackId = stepTrackId[s];
			minTime = signal_t;
		}
	}
	
	// If no step pass the threshold, getting the fastest signal - given by time + doca(s) / drift velocity
	if(trackId == -1)
	{
		for(unsigned int s=0; s<nsteps; s++)
		{
			G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z());
			signal_t = stepTime[s]/ns + DOCA.mag()/drift_velocity;
			if(signal_t < minTime)
			{
				trackId = stepTrackId[s];
				minTime = signal_t;
			}
		}
	}
	
	// Get Total Energy deposited of the fastest track
	//	double EtotF = 0;
	//	for(int s=0; s<nsteps; s++)
	//		if(stepTrackId[s] == trackId)
	//		{
	//			EtotF = EtotF + Edep[s];
	//		}
	

	// Distance 1: smear by 300 um
	// Drift Time 1: using 50 um/ns
	
	// Finding DOCA
	double doca    = 10000;
	double LR      = 0;
	
	// calculating LR
	for(unsigned int s=0; s<nsteps; s++)
	{
		G4ThreeVector DOCA(0, Lpos[s].y() + ylength - WIRE_Y, Lpos[s].z());
		if(DOCA.mag() <= doca && stepTrackId[s] == trackId )
		{
			doca = DOCA.mag();
			if(DOCA.y() >=0 ) LR = 1;
			else  LR = -1;
			
		}
	}
	
	double sdoca  = -1;
	sdoca = fabs(CLHEP::RandGauss::shoot(doca, 0.3));  ///< smeared by 300 microns for now
	double time1  = doca/drift_velocity;
	double stime1 = sdoca/drift_velocity;
	
	dgtz["hitn"]       = hitn;
	dgtz["sector"]     = sector;
	dgtz["superlayer"] = SL;
	dgtz["layer"]      = Layer;
	dgtz["wire"]       = nwire;
	dgtz["LR"]         = LR;
	dgtz["doca"]       = doca;
	dgtz["sdoca"]      = sdoca;
	dgtz["time"]       = time1;
	dgtz["stime"]      = stime1;
		
	return dgtz;
}

vector<identifier>  dc_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	double mini_stagger_shift;
	
	int SL                    = yid[1].id;
	int Layer                 = yid[2].id;
	
	mini_stagger_shift = 0;
	if(SL > 2 && SL <= 4)     // shift only for region 2
	{
		mini_stagger_shift = mini_stagger_shift_R2;
		if(Layer == 2 || Layer == 4 || Layer == 6)  mini_stagger_shift *= -1;
	}
	else if(SL > 4)           // shift only for region 3
	{
		mini_stagger_shift = mini_stagger_shift_R3;
		if(Layer == 2 || Layer == 4 || Layer == 6)  mini_stagger_shift *= -1;
	}
	
	
	if(mini_stagger_shift)
	cout << mini_stagger_shift << " mini stagger " << endl;
	
	G4StepPoint   *prestep   = aStep->GetPreStepPoint();
	G4StepPoint   *poststep  = aStep->GetPostStepPoint();
	G4VTouchable* TH = (G4VTouchable*) aStep->GetPreStepPoint()->GetTouchable();
	
	string         name    = TH->GetVolume(0)->GetName();                                    ///< Volume name
	G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction
	G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);
	
	
	double ylength = Detector.dimensions[3];  ///< G4Trap Semilength
	double deltay  = 2.0*ylength/NWIRES;
	double loc_y   = Lxyz.y() + ylength - mini_stagger_shift;      ///< Distance from bottom of G4Trap - modified by ministaggger
	
	int nwire = (int) floor(loc_y/deltay);
	if(nwire <= 0 )  nwire = 1;
	if(nwire == 113) nwire = 112;
	
	yid[3].id = nwire;
	
	if(fabs( (nwire+1)*deltay - loc_y ) < fabs( nwire*deltay - loc_y ) && nwire != 112 )
		yid[3].id = nwire + 1;
	
	/*	if(nwire > NWIRES) cout << " SuperLayer: " << SL << "      layer: "        << Layer
	 << "       wire: " << nwire     << "      length: "       << ylength
	 << "         ly: " << Lxyz.y()  << "      deltay: "       << deltay
	 << "      loc_y: " << loc_y     << "      nwire*deltay: " << fabs( nwire*deltay - loc_y ) << endl;*/
	
	yid[3].id_sharing = 1;
	
	return yid;
}



map< string, vector <int> >  dc_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}







