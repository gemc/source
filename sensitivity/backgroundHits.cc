// gemc
#include "backgroundHits.h"

// c++
#include <fstream>

// mlibrary
#include "gstring.h"
using namespace gstring;

// initialize from string
BackgroundHit::BackgroundHit(vector<string> hitsData)
{
	
}





// initialize map from filename
GBackgroundHits::GBackgroundHits(string filename, int verbosity)
{

	backgroundHitMap = new map<string, vector<BackgroundHit*> >;

	ifstream bgif(filename.c_str());

	if(bgif.good()) {

	} else {
		cout << " Warning: background file " << filename << " could not be opened. No background events will be merged." << endl;
		return;
	}

	if(verbosity > 0) cout << " Loading background hits from " << filename << endl;

	// file is good, loading hits
	while(!bgif.eof())
	{
		string bgline;
		getline(bgif, bgline);

		if(!bgline.size())
			continue;

		vector<string> hitsData = getStringVectorFromString(bgline);

		string systemEventNumber = hitsData[0] + "____" + hitsData[1];

		// load hits from string
		BackgroundHit *thisHit = new BackgroundHit(hitsData);

		// hits
		if(backgroundHitMap->find(systemEventNumber) == backgroundHitMap->end()) {

		}

	}




}






