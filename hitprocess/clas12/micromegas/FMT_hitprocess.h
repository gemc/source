#ifndef FMT_HITPROCESS_H
#define FMT_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"
#include "fmt_strip.h"


// Class definition

/// \class FMT_HitProcess
/// <b> FMT_HitProcess </b>\n\n
/// The forward micromegas hit process routine includes:\n
/// - conversion of the deposited energy into a number of electrons: the deposited energy provided by Geant4 is
/// used to estimate the average number of electrons created in the MM sensitive volume (knowing the mean
/// ionization potential of the gas). If the deposited energy is smaller than the ionization potential
/// (roughly 25 eV), then the hit is discarded
/// - generation of a transverse diffusion: for each electron created at the previous step, a transverse
/// diffusion is generated on a Gaussian. The width of this Gaussian is obtained by the standard formula for
/// diffusion, i.e. proportional to the square root of the distance between the interaction point and the micro-mesh.
/// Because of the focusing magnetic field, the amplitude of the transverse diffusion is considerably reduced
/// compared to the BMT.
/// - determination of the strip: the final transverse position of each electron (after transverse diffusion)
/// is converted into a strip number (the closest one).
/// - determination of the energy collected on the closest strip: if N electrons have been created for the
/// deposited energy Edep, the energy associated to the current strip is Edep/N. Note that this assumes no
/// gain fluctuation. A more realistic treatment (not yet implemented) would be to generate a gain for each
/// electron, using a Polya distribution.
/// - check of the acceptance: if the strip number is negative or larger than the total number of strips,
/// it is discarded (associated value is -1).
/// \author \n S. Procureur, G. Charles, M. Ungaro
/// \author mail: sebastien.procureur@cea.fr, gabriel.charles@cea.fr, ungaro@jlab.org\n\n\n

class FMT_HitProcess : public HitProcess
{
public:

	~FMT_HitProcess(){;}

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
	static HitProcess *createHitClass() {return new FMT_HitProcess;}

private:

	// constants initialized with initWithRunNumber
	static fmtConstants fmtc;

	void initWithRunNumber(int runno);
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

};

#endif
