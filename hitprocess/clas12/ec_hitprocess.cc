// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "ec_hitprocess.h"

static ecConstants initializeECConstants(int runno)
{
	ecConstants ecc;
	ecc.runNo = 0;
	
	ecc.NSTRIPS             = 36;
	ecc.attlen              = 3760.; // Attenuation Length (mm)
	ecc.TDC_time_to_channel = 20.;   // conversion from time (ns) to TDC channels.
	ecc.ECfactor            = 3.5;   // number of p.e. divided by the energy deposited in MeV; value taken from gsim. see EC NIM paper table 1.
	ecc.TDC_MAX             = 4095;  // max value for EC tdc.
	ecc.ec_MeV_to_channel   = 10.;   // conversion from energy (MeV) to ADC channels
	
	return ecc;
}

void ec_HitProcess::initWithRunNumber(int runno)
{
	if(ecc.runNo != runno)
	{
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		ecc = initializeECConstants(runno);
		ecc.runNo = runno;
	}
}


// Process the ID and hit for the EC using EC scintillator slab geometry instead of individual strips.
map<string, double> ec_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// get sector, stack (inner or outer), view (U, V, W), and strip.
	int sector = identity[0].id;
	int stack  = identity[1].id;
	int view   = identity[2].id;
	int strip  = identity[3].id;
	trueInfos tInfos(aHit);

	
	// Get scintillator mother volume dimensions (mm)
	double pDy1 = aHit->GetDetector().dimensions[3];  ///< G4Trap Semilength.
	double pDx2 = aHit->GetDetector().dimensions[5];  ///< G4Trap Semilength.
	double BA   = sqrt(4*pow(pDy1,2) + pow(pDx2,2)) ;
	
	vector<G4ThreeVector> pos  = aHit->GetPos();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	
	// Get Total Energy deposited
	double Etota = 0;
	double latt = 0;
	vector<G4double> Edep = aHit->GetEdep();
	
	for(unsigned int s=0; s<tInfos.nsteps; s++)
	{
		if(ecc.attlen>0)
		{
			double xlocal = Lpos[s].x();
			double ylocal = Lpos[s].y();
			if(view==1) latt = xlocal+(pDx2/(2.*pDy1))*(ylocal+pDy1);
			if(view==2) latt = BA*(pDy1-ylocal)/2./pDy1;
			if(view==3) latt = BA*(ylocal+pDy1-xlocal*2*pDy1/pDx2)/4/pDy1;
			Etota = Etota + Edep[s]*exp(-latt/ecc.attlen);
		}
		else
		{
			Etota = Etota + Edep[s];
		}
	}
	
	
	// initialize ADC and TDC
	int ADC = 0;
	int TDC = ecc.TDC_MAX;
	
	// simulate the adc value.
	if (Etota > 0)
	{
		double EC_npe = G4Poisson(Etota*ecc.ECfactor); //number of photoelectrons
		//  Fluctuations in PMT gain distributed using Gaussian with
		//  sigma SNR = sqrt(ngamma)/sqrt(del/del-1) del = dynode gain = 3 (From RCA PMT Handbook) p. 169)
		//  algorithm, values, and comment above taken from gsim.
		double sigma = sqrt(EC_npe)/1.22;
		double EC_MeV = G4RandGauss::shoot(EC_npe,sigma)*ecc.ec_MeV_to_channel/ecc.ECfactor;
		if (EC_MeV <= 0) EC_MeV = 0.0; // guard against weird, rare events.
		ADC = (int) EC_MeV;
	}
	
	// simulate the tdc.
	TDC = (int) (tInfos.time*ecc.TDC_time_to_channel);
	if (TDC > ecc.TDC_MAX) TDC = ecc.TDC_MAX;
	
	dgtz["hitn"]   = hitn;
	dgtz["sector"] = sector;
	dgtz["stack"]  = stack;
	dgtz["view"]   = view;
	dgtz["strip"]  = strip;
	dgtz["ADC"]    = ADC;
	dgtz["TDC"]    = TDC;
		
	return dgtz;
}

vector<identifier>  ec_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	
	int view        = yid[2].id; // get the view: U->1, V->2, W->3
	
	G4StepPoint   *prestep   = aStep->GetPreStepPoint();
	G4StepPoint   *poststep  = aStep->GetPostStepPoint();
	
	G4ThreeVector   xyz    = poststep->GetPosition();                                        ///< Global Coordinates of interaction
	G4ThreeVector  Lxyz    = prestep->GetTouchableHandle()->GetHistory()                     ///< Local Coordinates of interaction
	->GetTopTransform().TransformPoint(xyz);
	
	double xlocal = Lxyz.x();
	double ylocal = Lxyz.y();
	
	double pDy1 = Detector.dimensions[3];  ///< G4Trap Semilength. these points are needed to get the strip number.
	double pDx2 = Detector.dimensions[5];  ///< G4Trap Semilength.
	
	// Some EC geometry:
	//
	//    B ---------------------- C
	//       \                  /
	//        \                /
	//         \              /
	//          \            /            face of EC sector looking from the target; z-axis is out
	//           \     x    /
	//            \        /                      | y
	//             \      /                       |       coordinates: origin is at the
	//              \    /                 _______|       geometric center of the triangle.
	//               \  /                  x
	//                \/
	//                 A
	//
	//   U - strips parallel to BC, start numbering at A
	//   V - strips parallel to AB, start numbering at C
	//   W - strips parallel tp CA, start numbering at B
	//
	// other points: D - point where line from C crosses AB at right angle.
	//               F - CF is the component of the hit along CD.
	//               G - point where line from B crosses AC at right angle.
	//               H - BH is the component of the hit position along BG.
	//
	
	// get the strip.
	double BA = sqrt(4*pow(pDy1,2) + pow(pDx2,2)) ;
	//double lu=xlocal+(pDx2/(2.*pDy1))*(ylocal+pDy1);
	//double lv=BA*(pDy1-ylocal)/2./pDy1;
	//double lw=BA*(ylocal+pDy1-xlocal*2*pDy1/pDx2)/4/pDy1;
	//cout<<"pDx2="<<pDx2<<" pDy1="<<pDy1<<endl;
	
	int strip = 0;
	if (view == 1)
	{
		strip = (int) floor((ylocal + pDy1)*ecc.NSTRIPS/(2*pDy1)) + 1;  // U view strip.
		//cout << "view=" << view << " strip=" << strip << " xloc=" << xlocal << " yloc=" << ylocal << " pDy1: " << pDy1 << " floor:" << (ylocal + pDy1)*ecc.NSTRIPS/(2*pDy1) << " NS: " << ecc.NSTRIPS << endl;
	}
	else if (view == 2)
	{
		double BHtop = (xlocal + pDx2)*(-2*pDy1) + (ylocal - pDy1)*(pDx2);
		double BH = abs(BHtop/BA); // component of vector from the apex of the W view to the hit along a line perpendicular to the strips.
		
		double BGtop =  4*pDx2*pDy1;
		double BG = abs(BGtop/BA); // maximum length of W view along a line perpendicular to the strips (height of the W view triangle).
		
		strip = (int) floor(BH*ecc.NSTRIPS/BG)+1 ; // W view strip.
		//cout<<"view="<<view<<" strip="<<strip<<" xloc="<<xlocal<<" yloc="<<ylocal<<endl;
	}
	else if (view == 3)
	{
		double CFtop = (xlocal - pDx2)*2*pDy1 + (ylocal - pDy1)*pDx2;
		double CF = abs(CFtop/BA);  // component of vector from the apex of the V view to the hit along a line perpendicular to the strips.
		
		double CDtop =  -4*pDx2*pDy1;
		double CD = abs(CDtop/BA);  // maximum length of V view along a line perpendicular to the strips (height of the V view triangle).
		
		strip = (int) floor(CF*ecc.NSTRIPS/CD)+1 ; // V view strip.
		//cout<<"view="<<view<<" strip="<<strip<<" xloc="<<xlocal<<" yloc="<<ylocal<<endl;
	}
	else
	{
		cout << "   EC Hit Process WARNING: No view found." << endl;
	}
	
	
	if (strip <=0 ) strip = 1;
	if (strip >=36) strip = 36;
	
	yid[3].id = strip;
	yid[3].id_sharing = 1;
	
	return yid;
}



map< string, vector <int> >  ec_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}



// this static function will be loaded first thing by the executable
ecConstants ec_HitProcess::ecc = initializeECConstants(1);











