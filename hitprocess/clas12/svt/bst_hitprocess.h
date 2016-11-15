#ifndef BST_HITPROCESS_H
#define BST_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition



// constants to be used in the digitization routine
class bstConstants
{
public:

	// database
	int    runNo;
	string variation;
	string date;
	string connection;

	// translation table
	TranslationTable TT;

	// add constants here

};




/// \class bst_HitProcess
/// <b> bst_HitProcess </b>\n\n
/// The barrel silicon tracker hit process routine include:\n
/// - strip determination
/// - charge sharing
/// \author \n M. Ungaro
/// \author mail: ungaro@jlab.org\n\n\n


class bst_HitProcess : public HitProcess
{
public:
	~bst_HitProcess(){;}

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
	static HitProcess *createHitClass() {return new bst_HitProcess;}

private:

	// constants initialized with initWithRunNumber
	static bstConstants bstc;

	void initWithRunNumber(int runno);
	
	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();
};

#endif
