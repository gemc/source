#include "backgroundHits.h"


// initialize from string
BackgroundHit::BackgroundHit(string init)
{
	
}





// initialize map from filename
GBackgroundHits::GBackgroundHits(string filename)
{
	ifstream bgif(filename.c_str());

	if(bgif.good()) {

	} else {
		cout << " Warning: background file " << filename << " could not be opened. No background events will be merged." << endl;
	}


	// file is good, loading hits
	while(!bgif.eof())
	{
		string bgline;
		getline(IN, bgline);

		if(!bgline.size())
			continue;

	}




}






