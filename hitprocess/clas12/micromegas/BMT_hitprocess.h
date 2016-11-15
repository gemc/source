#ifndef BMT_HITPROCESS_H
#define BMT_HITPROCESS_H 1

/// \file BMT_hitprocess.h
/// BMT_hitprocess.h defines the barrel micromegas hit process routine class.\n
/// \author \n S. Procureur, G. Charles, M. Ungaro
/// \author mail: sebastien.procureur@cea.fr, gabriel.charles@cea.fr, ungaro@jlab.org\n\n\n


// gemc headers
#include "HitProcess.h"
#include "bmt_strip.h"


// Class definition

/// \class BMT_HitProcess
/// <b> BMT_HitProcess </b>\n\n
/// The barrel micromegas hit process routine includes:\n
/// - determination of the sector: the tile number is calculated using the phi coordinate of the hit
/// - conversion of the deposited energy into a number of electrons: the deposited energy provided by Geant4
/// is used to estimate the average number of electrons created in the MM sensitive volume (knowing the mean
/// ionization potential of the gas). If the deposited energy is smaller than the ionization potential
/// (roughly 25 eV), then the hit is discarded
/// - generation of a transverse diffusion: for each electron created at the previous step, a transverse
/// diffusion is generated on a Gaussian. The width of this Gaussian is obtained by the standard formula for
/// diffusion, i.e. proportional to the square root of the distance between the interaction point and the micro-mesh.
/// - shift of the position due to the Lorentz angle: for "z" detectors, the presence of the longitudinal
/// magnetic field induces a Lorentz angle which shifts the position of the signal. This shift is calculated for
/// each electron, and added to the transverse diffusion estimated at the previous step. Here again,
/// the shift due to the Lorentz angle depends on the distance between the interaction point and the micro-mesh.
/// The Lorentz angle has been estimated to 20 degrees using a Garfield simulation.
/// - determination of the strip: the final transverse position of each electron (after transverse diffusion
/// and potential Lorentz angle) is converted into a strip number (the closest one).
/// - determination of the energy collected on the closest strip: if N electrons have been created for the
/// deposited energy Edep, the energy associated to the current strip is Edep/N. Note that this assumes no gain
/// fluctuation. A more realistic treatment (not yet implemented) would be to generate a gain for each electron,
/// using a Polya distribution.
/// - check of the acceptance: if the strip number is negative or larger than the total number of strips,
/// it is discarded (associated value is -1).
/// \author \n V. Ziegler, S. Procureur, G. Charles, M. Ungaro
/// \author mail: ziegler@jlab.org sebastien.procureur@cea.fr, gabriel.charles@cea.fr, ungaro@jlab.org\n\n\n

// Veronique Ziegler (Dec. 3 2015)
// M. Ungaro (Jan. 26 2016)

class BMT_HitProcess : public HitProcess
{
public:

	~BMT_HitProcess(){;}

	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);

	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);

	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);

	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);

	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);

	// creates the HitProcess
	static HitProcess *createHitClass() {return new BMT_HitProcess;}


private:

	// constants initialized with initWithRunNumber
	static bmtConstants bmtc;

	double fieldScale;

	void initWithRunNumber(int runno);

	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
};

#endif
