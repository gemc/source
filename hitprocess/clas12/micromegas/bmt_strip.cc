// gemc headers
#include "bmt_strip.h"
#include "Randomize.hh"

// c++ headers
#include <iostream>
#include <cmath>
#include <cstdlib>


// Veronique Ziegler (Dec. 3 2015)
// Note: this method only contains the constants for the third micromegas region,
// the constants for the other 2 are missing...
// M. Ungaro (Jan 26 2016)
// Updated to read constants at the beginning of each run

// the routine to find the strip
vector<double> bmt_strip::FindStrip(int layer, int sector, G4ThreeVector xyz, double Edep, bmtConstants bmtc)
{
	double x = xyz.x()/mm;
	double y = xyz.y()/mm;
	double z = xyz.z()/mm;

	double w_i = 25; //ionization potential assumed to be 25 eV
	int Nel = (int) (1e6*Edep/w_i); // the number of electrons produced

	// the return vector is always in pairs the first index is the strip number, the second is the Edep on the strip
	vector<double> strip_id;

	int strip = -1;
	if(Nel>0) // if the track deposited energy is greater than the assumed ionization potential digitize
	{
		for(int iel=0;iel<Nel;iel++)  // at this stage the same algorithm requiring looping over all e-'s is kept (slow...)
		{
			int num_region = (int) (layer+1)/2 - 1; 	// region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6;
			// start the loop and get the strip number from the layer, z and phi hit information
			double sigma =0;
			if(layer%2==0)
			{// C layer
				double z_i = bmtc.CRZZMIN[num_region] + bmtc.CRZOFFSET[num_region]; 						 // fiducial z-profile lower limit
				double z_f = bmtc.CRZZMIN[num_region] + bmtc.CRZOFFSET[num_region] + bmtc.CRZLENGTH[num_region];  // fiducial z-profile upper limit

				if(z>=z_i && z<=z_f)
				{
					sigma = getSigmaLongit(layer, x, y, bmtc);  //  longitudinal diffusion

					double z_d = (double) (G4RandGauss::shoot(z,sigma));	 //  shower profile in z generated using a gaussian distribution with width=sigma

					strip = getCStrip(layer, z_d, bmtc);
				}
			}

			if(layer%2==1)
			{// Z layer
				sigma = getSigmaAzimuth(layer, x, y, bmtc); //   diffusion  taking into account the Lorentz angle

				double Delta_rad = sqrt(x*x+y*y) - bmtc.CRZRADIUS[num_region] + bmtc.hStrip2Det;

				double phi = atan2(y,x);
				
				//shower profile in phi generated
				double phi_d = phi + ((G4RandGauss::shoot(0,sigma))/cos(bmtc.ThetaL) - Delta_rad*tan(bmtc.ThetaL))/bmtc.CRZRADIUS[num_region];

				strip = getZStrip(layer, phi_d, bmtc);
			}
			// if the strip is found (i.e. the hit is within acceptance
			if(strip != -1)
			{	// search the existing strip array to see if this strip is already stored, if yes, add contributing energy
				for(int istrip=0;istrip< (int) (strip_id.size()/2);istrip++)
				{
					if(strip_id[2*istrip]==strip)
					{// already hit strip - add Edep
						strip_id[2*istrip+1] = strip_id[2*istrip+1]+1./((double) Nel); // no gain fluctuation ;
						strip = -1; // not to use it anymore
					}
				}
				if(strip > -1) { // new strip, add it
					strip_id.push_back(strip);
					strip_id.push_back(1./((double) Nel)); // no gain fluctuation yet
				}
			}
			else
			{// not in the acceptance
				strip_id.push_back(-1);
				strip_id.push_back(1);
			}
		}
	}
	else
    { // Nel=0, consider the Edep is 0
        strip_id.push_back(-1);
        strip_id.push_back(1);
    }

    return strip_id;
}


//
// param layer
// param x x-coordinate of the hit in the lab frame
// param y y-coordinate of the hit in the lab frame
// return the sigma along the beam direction (longitudinal)
//
double bmt_strip::getSigmaLongit(int layer, double x, double y, bmtConstants bmtc)
{
	
	// sigma for C-detector
	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	return bmtc.SigmaDrift*sqrt((sqrt(x*x+y*y) - bmtc.CRCRADIUS[num_region] + bmtc.hStrip2Det)/bmtc.hDrift);
}

//
// param layer
// param x x-coordinate of the hit in the lab frame
// param y y-coordinate of the hit in the lab frame
// return the sigma along in the azimuth direction taking the Lorentz angle into account
//
double bmt_strip::getSigmaAzimuth(int layer, double x, double y,  bmtConstants bmtc)
{ // sigma for Z-detectors

	int num_region = (int) (layer+1)/2 - 1; ///< region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	double sigma = bmtc.SigmaDrift*sqrt((sqrt(x*x+y*y) - bmtc.CRZRADIUS[num_region] + bmtc.hStrip2Det)/bmtc.hDrift/cos(bmtc.ThetaL));

	return sigma;

}

//
// param layer the layer 1...6
// param angle the position angle of the hit in the Z detector
// return the Z strip as a function of azimuthal angle
//
int bmt_strip::getZStrip(int layer, double angle, bmtConstants bmtc)
{
	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	int num_detector =isInSector( layer,  angle, bmtc) - 1;
	if(num_detector==-1)
		return -1;

	if(angle<0)
		angle+=2*pi; // from 0 to 2Pi

	if(num_detector==1)
	{
		double angle_f=bmtc.CRCEDGE1[num_region][1] + (bmtc.CRCXPOS[num_region] + bmtc.CRCLENGTH[num_region])/bmtc.CRCRADIUS[num_region] - 2*pi;
		if(angle>=0 && angle<=angle_f)
			angle+=2*pi;
	}
	double strip_calc = ( (angle - bmtc.CRZEDGE1[num_region][num_detector])*bmtc.CRZRADIUS[num_region]-bmtc.CRZXPOS[num_region] - bmtc.CRZWIDTH[num_region]/2.)/(bmtc.CRZWIDTH[num_region]
																																																					 + bmtc.CRZSPACING[num_region]);
	
	strip_calc = (int)(round(strip_calc) );
	int strip_num = (int) floor(strip_calc);

	int value = strip_num + 1;

	if(value < 1 || value > bmtc.CRZNSTRIPS[num_region])
		value = -1;

	return value;
}
//
// param sector
// param layer
// param trk_z the track z position of intersection with the C-detector
// return the C-strip
//
int bmt_strip::getCStrip(int layer, double trk_z, bmtConstants bmtc)
{

	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	int strip_group = 0;
	int ClosestStrip =-1;
	// get group
	int len = bmtc.CRCGROUP[num_region].size();
	double Z_lowBound[len];
	double Z_uppBound[len];
	int NStrips[len];

	double zi= bmtc.CRCZMIN[num_region] + bmtc.CRCOFFSET[num_region];
	double z = trk_z - zi;

	Z_lowBound[0] = bmtc.CRCWIDTH[num_region][0]/2.; // the lower bound is the zMin+theOffset with half the width
	Z_uppBound[0] = Z_lowBound[0]
					   + (bmtc.CRCGROUP[num_region][0]-1)*(bmtc.CRCWIDTH[num_region][0]+ bmtc.CRCSPACING[num_region]);
	NStrips[0] = bmtc.CRCGROUP[num_region][0];
	for(int i =1; i< len; i++)
	{
		Z_lowBound[i] = Z_uppBound[i-1] + bmtc.CRCWIDTH[num_region][i-1]/2. + bmtc.CRCSPACING[num_region] + bmtc.CRCWIDTH[num_region][i]/2.;
		Z_uppBound[i] = Z_lowBound[i]   + (bmtc.CRCGROUP[num_region][i]-1)*(bmtc.CRCWIDTH[num_region][i] + bmtc.CRCSPACING[num_region]);

		NStrips[i] = NStrips[i-1] + bmtc.CRCGROUP[num_region][i];

		if(z>=Z_lowBound[i] && z<=Z_uppBound[i]) {
			strip_group = i;
			ClosestStrip = 1 + (int) (round(((z-Z_lowBound[strip_group])/(bmtc.CRCWIDTH[num_region][strip_group] + bmtc.CRCSPACING[num_region]))))+NStrips[i-1];

			len =i;
		}
	}
	return ClosestStrip;
}



// param layer the hit layer
// param strip the hit strip
// return the z position in mm for the C-detectors
//
// WARNING: This routine is not used?
//
double bmt_strip::CRCStrip_GetZ(int layer, int strip, bmtConstants bmtc)
{
	int num_strip = strip - 1;     			// index of the strip (starts at 0)
	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	//For CR6C, this function returns the Z position of the strip center
	int group=0;
	int limit = bmtc.CRCGROUP[num_region][group];
	double zc = bmtc.CRCZMIN[num_region] + bmtc.CRCOFFSET[num_region] + bmtc.CRCWIDTH[num_region][group]/2.;

	if (num_strip>0){
	  for (int j=1;j<num_strip+1;j++)
	  {
		zc += bmtc.CRCWIDTH[num_region][group]/2.;
		if (j>=limit)
		{ //test if we change the width
			group++;
			limit += bmtc.CRCGROUP[num_region][group];
		}
		zc += bmtc.CRCWIDTH[num_region][group]/2. + bmtc.CRCSPACING[num_region];
	  }
	}

	return zc; //in mm
}

// param sector the sector in CLAS12 1...3
// param layer the layer 1...6
// param strip the strip number (starts at 1)
// return the angle to localize the  center of strip
//
double bmt_strip::CRZStrip_GetPhi(int sector, int layer, int strip, bmtConstants bmtc)
{
	// Sector = num_detector + 1;
	// num_detector = 0 (region A), 1 (region B), 2, (region C)
	//For CRZ, this function returns the angle to localize the  center of strip "num_strip" for the "num_detector"
	int num_detector = getDetectorIndex(sector); 			// index of the detector (0...2): sector 1 corresponds to detector B, 2 to A,  3 to C in the geometry created by GEMC
	int num_strip = strip - 1;     							// index of the strip (starts at 0)
	int num_region = (int) (layer+1)/2 - 1; 				// region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	double angle = bmtc.CRZEDGE1[num_region][num_detector] + (bmtc.CRZXPOS[num_region]+(bmtc.CRZWIDTH[num_region]/2.
																													+ num_strip*(bmtc.CRZWIDTH[num_region]+bmtc.CRZSPACING[num_region])))/bmtc.CRZRADIUS[num_region];
	if (angle>2*pi) angle-=2*pi;
	return angle; //in rad
}

// not used (implemented for alternate algorithm not yet fully developed)
double bmt_strip::getEnergyFraction(double z0, double z, double sigma)
{
	double pdf_gaussian = (1./(sigma*sqrt(2*pi)))* exp( -0.5*((z-z0)/sigma)*((z-z0)/sigma) );
	return pdf_gaussian;
}

int bmt_strip::getDetectorIndex(int sector)
{
	//sector 1 corresponds to detector B, 2 to A,  3 to C
	// A is detIdx 0, B 1, C 2.
	int DetIdx = -1;
	if(sector == 1)
		DetIdx = 1;
	if(sector == 2)
		DetIdx = 0;
	if(sector == 3)
		DetIdx = 2;

	return DetIdx;
}


// not used yet (implemented for alternate algorithm not yet fully developed)
int bmt_strip::isInSector(int layer, double angle, bmtConstants bmtc)
{

	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	if(angle<0)
		angle+=2*pi; // from 0 to 2Pi
	double angle_pr = angle + 2*pi;

	double angle_i = 0; // first angular boundary init
	double angle_f = 0; // second angular boundary for detector A, B, or C init
	int num_detector = -1;

	for(int i = 0; i<3; i++) {

		angle_i = bmtc.CRCEDGE1[num_region][i] + bmtc.CRCXPOS[num_region] /bmtc.CRCRADIUS[num_region];
		angle_f = bmtc.CRCEDGE1[num_region][i] + (bmtc.CRCXPOS[num_region] + bmtc.CRCLENGTH[num_region])/bmtc.CRCRADIUS[num_region];

		if( (angle>=angle_i && angle<=angle_f) || (angle_pr>=angle_i && angle_pr<=angle_f) )
			num_detector=i;
	}
	return num_detector + 1;
}





