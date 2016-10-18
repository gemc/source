#include <vector>
using namespace std;

// geant4 includes
#include "G4ThreeVector.hh"

class bst_strip
{
public:
	double alpha;
	double pitch;

	double intersensors;       // gap between sensors in the same module, different layers
	vector<int>    NSensors;   // number of sensors by sector for each layer

	double DZ_inLength;   // size of the band of dead zones all around in the length of the card
	double DZ_inWidth;    // size of the band of dead zones all around in the width of the card
	double SensorLength;  // length of 1 Sensor (including dead area)
	double SensorWidth ;  // width 1 Sensor (including dead area)
	double StripStart;    // the first sensitive strip starts at stripStart + 1/2 of the pitch
	int Nstrips;          // Number of strips for 1 card (New Design)

	void fill_infos();

	vector<double> FindStrip( int layer, int sector, int isens, G4ThreeVector Lxyz);   // Strip Finding Routine

};
