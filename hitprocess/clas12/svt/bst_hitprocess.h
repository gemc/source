#ifndef BST_HITPROCESS_H
#define BST_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition

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
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
		
		// creates the HitProcess
		static HitProcess *createHitClass() {return new bst_HitProcess;}
};

#endif
