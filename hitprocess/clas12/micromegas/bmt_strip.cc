// gemc headers
#include "bmt_strip.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

void bmt_strip::fill_infos()
{
 // all dimensions are in mm
    Pi  = 3.14159265358;

    interlayer     = 15.0;
    Nsector        = 3;
    hDrift         = 3.0;
    hStrip2Det     = hDrift/2.;
    sigma_td_max   = 0.4;
    theta_L        = Pi*20./180.;
    w_i            = 25.0;


    interStripZ = 0.20;
    interStripC = 0.16;

 // pitch CR4Z -> CR4 = central tracker region 4
    pitchZ4 = 0.487; // and interpitch=0.201
 
 // pitch CR5Z -> CR5 = central tracker region 5
    pitchZ5 = 0.49; // none existing value. picked a rabdom value compatible with the geometry
 
 // pitch CR6Z -> CR6 = central tracker region 6
    double strip_widthZ = 0.328;
    pitchZ6 = interStripZ + strip_widthZ;
 
 // pitch CRC --> for the C-detectors, the strips are in bunches of equal pitches
    double CR4C_width[13]={0.345,0.28,0.225,0.175,0.17,0.21,0.26,0.31,0.37,0.44,0.515,0.605,0.7};     // width of the corresponding group of strips
    double CR4C_group[13]={32,32,32,32,624,32,32,32,32,32,32,32,896};                              	  // the number of strips with equal pitches
    // For CR5, no existing value. picked a random value compatible with the geometry
    double CR5C_width[1]={0.253};                                                                     // width of the corresponding group of strips
    double CR5C_group[1]={1024};                                                                      // the number of strips with equal pitches
    double CR6C_width[14]={0.38,0.32,0.27,0.23,0.17,0.18,0.22,0.25,0.29,0.33,0.37,0.41,0.46,0.51};    // width of the corresponding group of strips
    double CR6C_group[14]={32,32,32,32,704,64,32,32,32,32,32,32,32,32};                               // the number of strips with equal pitches
    int arraySize[] = {sizeof(CR4C_group) / sizeof(CR4C_group[0]),sizeof(CR5C_group) / sizeof(CR5C_group[0]),sizeof(CR6C_group) / sizeof(CR6C_group[0])}; // size of arrays describing the C-detectors for a given MM region

    // pitch CR4C for first region of macromegas (CRC4)
    for(int i =0; i<arraySize[0]; i++) {
    	widthC4.push_back(CR4C_width[i]);
        pitchC4.push_back(interStripC+CR4C_width[i]);
        nbunchC4.push_back(CR4C_group[i]);
    }
    // pitch CR5C for second region of macromegas (CRC5)
    for(int i =0; i<arraySize[1]; i++) {
    	widthC5.push_back(CR5C_width[i]);
        pitchC5.push_back(interStripC+CR5C_width[i]);
        nbunchC5.push_back(CR5C_group[i]);
    }
    // pitch CR6C for third region of macromegas (CRC6)
    for(int i =0; i<arraySize[2]; i++) {
    	widthC6.push_back(CR6C_width[i]);
        pitchC6.push_back(interStripC+CR6C_width[i]);
        nbunchC6.push_back(CR6C_group[i]);
    }
    // z of the upstream part of the layer
    Z0.push_back(-127.820); Z0.push_back(-127.820);
    Z0.push_back(-148.740); Z0.push_back(-148.740);
    Z0.push_back(-169.710); Z0.push_back(-169.710);
 
    // total z length of the layer of the layer
    DZ.push_back(372.75); DZ.push_back(372.75);
    DZ.push_back(423.99); DZ.push_back(423.99);
    DZ.push_back(444.96); DZ.push_back(444.96);

    // radii of layers
    R.push_back(145.731); R.push_back(R[0]+interlayer);
    R.push_back(175.731); R.push_back(R[2]+interlayer);
    R.push_back(205.731); R.push_back(R[4]+interlayer);

    // Number of strips (depends on radius, dead zones, and pitch!)
    // layers 0,2,4 are Z detectors, i.e. measure phi (strips in z direction)
    // layers 1,3,5 are C detectors, i.e. measure z (phi strips at z)
    Nstrips.push_back(14*64);         // CR4Z
    Nstrips.push_back(10*64);         // CR4C
    Nstrips.push_back(11*64);         // CR5Z
    Nstrips.push_back(16*64);        // CR5C
    Nstrips.push_back(768);          // CR6Z
    Nstrips.push_back(1152);         // CR6C

 
    // dead angles
    Inactivtheta.push_back((20/R[0])*(180./Pi));
    Inactivtheta.push_back((20/R[1])*(180./Pi));
    Inactivtheta.push_back((20/R[2])*(180./Pi));
    Inactivtheta.push_back((20/R[3])*(180./Pi));
    Inactivtheta.push_back((20/R[4])*(180./Pi));
    Inactivtheta.push_back((20/R[5])*(180./Pi));
 
    // dead zones
    DZ_inLength = 0; //for CRnZ
    DZ_inWidth  = 0; //for CRnC
    DZ4_inLength = ((R[1]*(120-Inactivtheta[1])*Pi/180)-(Nstrips[1]*pitchZ4+interStripZ))/2; //for CR4Z (((160.731*(120-7.129)*Pi/180)-(10*64*0.487)+0.201))/2 = 2.377
    DZ4_inWidth  = (DZ[0]-372.48-interStripC)/2; //for CR4C // (372.75 - 372.48 - 0.16)/2 = 0.055
    DZ5_inLength = ((R[2]*(120-Inactivtheta[2])*Pi/180)-(Nstrips[2]*pitchZ4+interStripZ))/2; //for CR5Z (((175.731*(120-6.521)*Pi/180)-(11*64*0.49)+0.201))/2 = 1.444
    DZ5_inWidth  = (DZ[3]-422.912-interStripC)/2; //for CR5C // (423.99 - 422.912 - 0.16)/2 = 0.459
    DZ6_inLength = ((R[4]*(120-Inactivtheta[4])*Pi/180)-(Nstrips[4]*pitchZ4+interStripZ))/2; //for CR6Z (((205.731*(120-5.57)*Pi/180)-(12*64*0.526)+0.201))/2 = 1.779
    DZ6_inWidth  = (DZ[5]-444.8-interStripC)/2; //for CR6C // (444.96 - 44.8 - 0.16)/2 = 0
 
    // mid angle of the sector
    MidTile.push_back(0);     MidTile.push_back(0);
    MidTile.push_back(0);     MidTile.push_back(0);
    MidTile.push_back(0);     MidTile.push_back(0);
}
/**
 * Method to get digi hits based on x,y,z position and E
 * layers from 0 to 5
 *   layers 0,2,4 are Z detectors (inner)
 *   layers 1,3,5 are C detectors (outer)
 */
vector<double>  bmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
        // the return vector is always in pairs.
        // The first number is the ID,
        // the second number is the sharing percentage
        vector<double> strip_id;
        int NbStrips =0 ;
        // dead zones
        if(layer == 0 || layer == 1){
                DZ_inLength = DZ4_inLength;
                DZ_inWidth  = DZ4_inWidth;
        }
        else if(layer == 2 || layer == 3){
                DZ_inLength = DZ5_inLength;
                DZ_inWidth  = DZ5_inWidth;
        }
        else if(layer == 4 || layer == 5){
                DZ_inLength = DZ6_inLength;
                DZ_inWidth  = DZ6_inWidth;
        }


        // 1st define phi of the mean hit point
        double phi;
        if(x>0 && y>=0) phi = atan(y/x);
        else if(x>0 && y<0) phi = 2.*Pi+atan(y/x);
        else if(x<0) phi = Pi+atan(y/x);
        else if(x==0 && y>0) phi = Pi/2.;
        else if(x==0 && y<0) phi = 3.*Pi/2.;
        else phi = 0; // x = y = 0, phi not defined

        // now find the tile number (transverse diffusion will not change that)
        int ti=0;
        double theta_tmp=0;
        for(int t=0; t<Nsector; t++)
        {
                theta_tmp = MidTile[layer]+2.*t*Pi/Nsector;
                if(theta_tmp>2.*Pi) theta_tmp = theta_tmp - 2.*Pi;
                if(theta_tmp<0) theta_tmp = theta_tmp + 2.*Pi;
                if(fabs(phi-theta_tmp)<Pi/Nsector || fabs(2.*Pi-fabs(phi-theta_tmp))<Pi/Nsector) ti=t; // gives tile #
        }

        double phiij = MidTile[layer]+2.*ti*Pi/Nsector;
        if(phi-phiij<=-Pi/Nsector) phi = phi+2.*Pi;
        if(phi-phiij>Pi/Nsector) phi = phi-2.*Pi;
        if(phi-phiij<=-Pi/Nsector || phi-phiij>Pi/Nsector) cout << "WARNING: incorrect phi value in BMT: " << phi*180./Pi << " vs " << phiij*180./Pi << endl;

        // now compute the sigma of the (transverse) dispersion for this interaction
        if(layer==1 || layer==3 || layer==5) sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/hDrift); // "C" det, transverse diffusion grows with square root of distance
        else sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/(cos(theta_L)*hDrift)); // same, but "Z" detectors, so Lorentz angle makes drift distance longer by 1./cos(theta_L) . Means sigma_td can be larger than sigma_td_max
        cout<<" x "<<x<<" y "<<y<<" z "<<z<<endl;
        if(Edep>0)
        {
			NbStrips = Nstrips[layer];
			cout<<" NbStrips "<<NbStrips<<endl;
			if(layer%2==1)
			{ //  for "C" layers, i.e. measuring z
				vector<double> pitchC;
				vector<double> widthC;
				vector<int>nbunchC;
				int NbStrips;

				if(layer==1) {
					pitchC = pitchC4;
					widthC = widthC4;
					nbunchC = nbunchC4;
				}
				if(layer==3) {
					pitchC = pitchC5;
					widthC = widthC5;
					nbunchC = nbunchC5;
				}
				if(layer==5) {
					pitchC = pitchC6;
					widthC = widthC6;
					nbunchC = nbunchC6;
				}

				double z_min = z-3*sigma_td; // minimum z in range
				double z_max = z+3*sigma_td; // maximum z in range
				cout<<" z_min "<<z_min <<" z_max "<<z_max<<endl;
				double lowerBound = Z0[layer]+DZ_inWidth;
				double upperBound = Z0[layer]+DZ[layer]-DZ_inWidth;

				if(z>z_min && z_min<lowerBound)  	// the z_min falls outside of fiducial area
					z_min=lowerBound;				// move z_min to the lower edge of the fiducial area
				if(z<z_max && z_max>upperBound)  	// the z_max falls outside of fiducial area
					z_max=upperBound;				// move z_max to the upper edge of the fiducial area

				if(z_min>=lowerBound && z_max<=upperBound)
				{
					int min_strip = getNearestCstrip(z_min, layer, pitchC, nbunchC, NbStrips, DZ_inWidth);
					int max_strip = getNearestCstrip(z_max, layer, pitchC, nbunchC, NbStrips, DZ_inWidth);
					cout<<" min_strip "<<min_strip <<" max_strip "<<max_strip<<endl;
					for(int s = min_strip; s<=max_strip; s++)
					{
						double Cz = getZasfcnCstrip(s, layer, pitchC, widthC, nbunchC);
						double f =  getEnergyFraction(z, Cz, sigma_td);
						cout<<" strip "<<s <<"  f "<< f<<endl;
						if(f>0.05)
						{// 5% of total Edep cut off
							strip_id.push_back(s);
							strip_id.push_back(Edep*f);
						}
					}
				} else {
					strip_id.push_back(-1); // not in fiducial
					strip_id.push_back(0);
				}

			}
			if(layer%2==0)
			{ //  for "Z" layers, i.e. measuring phi
				int pitchZ =0;

				if(layer==0)
					pitchZ = pitchZ4;
				if(layer==2)
					pitchZ = pitchZ5;
				if(layer==4)
					pitchZ = pitchZ6;
                cout<<" layer "<< layer <<" pitchZ "+pitchZ<<endl;
				double phi_min = phi + ((-3*sigma_td)/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer];
				double phi_max = phi + ((3*sigma_td)/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer];

				double lowerBound = phiij-Pi/Nsector+DZ_inLength/R[layer] ;
				double upperBound = phiij+Pi/Nsector-DZ_inLength/R[layer] ;

				if(phi>phi_min && phi_min<lowerBound)  	// the phi_min falls outside of fiducial area
					phi_min=lowerBound;					// move phi_min to the lower edge of the fiducial area
				if(phi<phi_max && phi_max>upperBound)  	// the phi_max falls outside of fiducial area
					phi_max=upperBound;					// move phi_max to the upper edge of the fiducial area

				if(phi_min>=lowerBound && phi_max<=upperBound)
				{
					int min_strip = getNearestZstrip(layer, phi_min, phiij, pitchZ, DZ_inLength);
					int max_strip = getNearestZstrip(layer, phi_max, phiij, pitchZ, DZ_inLength);
					cout<<" phi "<<phi <<" phi_min "<<phi_min<<" phi_max "<<phi_max<<endl;
					cout<<" min_strip "<<min_strip<<" max_strip "<<max_strip<<endl;
					for(int s = min_strip; s<=max_strip; s++)
					{
						double Cphi = getPhiasfcnCstrip(s, layer, phi_min, phiij, pitchZ, DZ_inLength);
						double f =  getEnergyFraction(0, Cphi, sigma_td);
						if(f>0.05)
						{// 5% of total Edep cut off
							cout<<" strip "<<s <<"  f "<< f<<endl;
							strip_id.push_back(s);
							strip_id.push_back(Edep*f);
						}
					}
				} else {
					strip_id.push_back(-1); // not in fiducial
					strip_id.push_back(0);
				}
			}
        }
		else
		{ // Nel=0, consider the Edep is 0
				strip_id.push_back(-1);
				strip_id.push_back(0);
		}

        return strip_id;
}
int bmt_strip::getNearestCstrip(double z, int layer, vector<double> pitchC, vector<int>nbunchC, int NbStrips, double DZ_inWidth){

	int arraySize = nbunchC.size();
	int ClosestStrip =-1;
	double lowerBound = Z0[layer]+DZ_inWidth;
	double upperBound = Z0[layer]+DZ[layer]-DZ_inWidth;
	if(z==lowerBound)
		ClosestStrip = 1;
	if(z==upperBound)
			ClosestStrip = NbStrips;
	// if this z is within the first group of equal-pitch strips then get its closest strip
	if(z-lowerBound>0 && z-lowerBound<nbunchC[0]*pitchC[0])
			ClosestStrip = (int) (floor(((z-lowerBound)/pitchC[0]))+0.5);
	// the closest strip is the effective z divided by the pitch +0.5
	// if there is only one group -- stop here...
	// otherwise loop over the number of groups of equal-pitch
	// strips and find the corresponding nearest strip from z
	int nstripsPrevGrp =0;
	if(arraySize>1)
		for(int i =1; i<arraySize; i++) {
			nstripsPrevGrp+=i*nbunchC[i-1];
			if(z-lowerBound>nbunchC[i-1]*pitchC[i-1] && z-lowerBound<nbunchC[i]*pitchC[i])
					ClosestStrip = (int) (floor(((z-lowerBound)/pitchC[i]))+0.5+nstripsPrevGrp);
		}
	return ClosestStrip;
}
/**
 * A method to return the z position of a given C strip...
 */
double bmt_strip::getZasfcnCstrip(int strip, int layer, vector<double> pitchC, vector<double> widthC, vector<int>nbunchC){

	int num_strip = strip - 1;     			// index of the strip (starts at 0)

	//For CRC, this function returns the Z position of the strip center
	int group=0;
	int limit=nbunchC[group];
	double zc=Z0[layer]+widthC[group]/2.;

	if (num_strip>0){
	  for (int j=1;j<num_strip+1;j++){
		zc+=widthC[group]/2.;
		if (j>=limit) { //test if we change the width
			group++;
			limit+=widthC[group];
		}
		zc+=widthC[group]/2.+interStripC;
	  }
	}
	return zc;


}
int bmt_strip::getNearestZstrip(int layer, double phi, double phiij, double pitchZ, double DZ_inLength) {
	cout<<"layer "<<layer<<" phi "<<phi<<" phiij "<<phiij<<" pitchZ "<< pitchZ << endl;
	return (int) (floor(((R[layer]/pitchZ)*(phi-phiij+Pi/Nsector - (Inactivtheta[layer]/2.)*Pi/180. - DZ_inLength/R[layer]))+0.5));
}

double bmt_strip::getPhiasfcnCstrip(int s, int layer, double phi, double phiij, double pitchZ, double DZ_inLength) {
	return (s - 0.5)*pitchZ/R[layer]+phiij-Pi/Nsector + (Inactivtheta[layer]/2.)*Pi/180. + DZ_inLength/R[layer];
}
double bmt_strip::getEnergyFraction(double z0, double z, double sigma){
	double pdf_gaussian = (1./(sigma*sqrt(2*Pi)))* exp( -0.5*((z-z0)/sigma)*((z-z0)/sigma) );
	return pdf_gaussian;
}
