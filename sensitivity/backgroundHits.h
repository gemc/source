#ifndef backgroundHit_H
#define backgroundHit_H 1

// gemc
#include "identifier.h"

// c++
#include <map>
using namespace std;

class BackgroundHit {

public:
	// constructor from file directly
	BackgroundHit(string init);

private:
	
	double energy; // in MeV
	double timeAtElectronics; //
	double npheD; // number of photoelectrons detected
	
	vector<identifier> identity;
};


class GBackgroundHits {

	// initialize map from file
public:
	GBackgroundHits() = default;
	GBackgroundHits(string filename, int verbosity = 0);

	vector<BackgroundHit*> getBackgroundForSystem(string system);

private:
	// the key is detectorName____eventNumber
	map<string, vector<BackgroundHit*> > *backgroundHitMap;
	
};


#endif
