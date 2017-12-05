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
	GBackgroundHits(string filename);
	
	
private:
	
	map<string, vector<BackgroundHit*> > *backgroundHitMap;
	
};


#endif
