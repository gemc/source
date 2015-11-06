#include <vector>
using namespace std;

class fmt_strip
{
	public:
		double pitch;
		double Pi;

		double interlayer;       // distance between 2 layers of a superlayer
		double intersuperlayer;  // distance between 2 superlayers
		double R_max;            // outer radius of strip part
		double R_min;            // inner radius of strip part

		double Z_1stlayer;       // z position of the 1st layer

		vector<double> Z0;       // z of the upstream part of the layer
		vector<double> alpha;    // strip angles of layers
		int N_str;               // Number of strips for 1 layer
		double hDrift;           // Size of the drift gap
		double hStrip2Det;       // distance between strips and the middle of the conversion gap (~half the drift gap)
		double sigma_td_max;     // maximum value of the transverse diffusion
		double sigma_td;         // current value of the transverse diffusion
		double w_i;              // mean ionization potential
		int Nel;                 // number of electrons (Nt) for a given hit
		double y_real; 			 // y position of the hit after transverse diffusion
		double x_real;			 // x position of the hit after transverse diffusion
		void fill_infos();

		vector<double> FindStrip( int layer, int sector, double x, double y, double z, double Edep);   // Strip Finding Routine

};
