// gemc headers
#include "bmt_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>
#include <cstdlib>

void bmt_strip::fill_infos()
{
 // all dimensions are in mm
    Pi  = 3.14159265358;
    SigmaDrift = 0.4;
    hDrift = 3.0;
    hStrip2Det = hDrift/2.;
    FieldFlag = 0;
    ThetaL = (double) FieldFlag*toRadians(20.);

    // DEFINING THE Z DETECTORS
    CRZRADIUS[2]	=	205.8;
	CRZNSTRIPS[2]	=	768;
	CRZSPACING[2]	=	0.2;
	CRZWIDTH[2]		=	0.328;
	CRZLENGTH[2]	=	444.88;
	CRZZMIN[2]		=	-421.75;
	CRZZMAX[2]		=	290.25;
	CRZOFFSET[2]	=	252.1;
    CRZEDGE1.resize(3); //3 regions
    for (int i = 0; i <3 ; ++i)
    	CRZEDGE1[i].resize(3);
    CRZEDGE2.resize(3); //3 regions
	for (int i = 0; i <3 ; ++i)
		CRZEDGE2[i].resize(3);

    double ZEdge1[] = {toRadians(30.56),toRadians(270.56),toRadians(150.56)};
    double ZEdge2[] = {toRadians(149.44),toRadians(29.44),toRadians(269.44)};

    // Assume the edge boundaries are the same for all regions until get final geometry
    for (int i = 0; i <3 ; ++i)
    {
    	for (int j = 0; j <3 ; ++j)
    	{
			CRZEDGE1[i][j] = ZEdge1[j];
			CRZEDGE2[i][j] = ZEdge2[j];
		}
    }

    CRZXPOS[2]		=	10.547;
    // DEFINING THE C DETECTORS
	CRCRADIUS[2]	=	220.8;
	CRCNSTRIPS[2]	=	1152;
	CRCLENGTH[2]	=	438.6;
	CRCSPACING[2]	=	0.16;
	CRCZMIN[2]		=	-421.75;
	CRCZMAX[2]		=	290.25;
	CRCOFFSET[2]	=	252.18;

    // pitch CRC --> for the C-detectors, the strips are in bunches of equal pitches
    double CR4C_width[13]={0.345,0.28,0.225,0.175,0.17,0.21,0.26,0.31,0.37,0.44,0.515,0.605,0.7};     // width of the corresponding group of strips
    double CR4C_group[13]={32,32,32,32,624,32,32,32,32,32,32,32,896};                              	  // the number of strips with equal pitches
    // For CR5, no existing value. picked a random value compatible with the geometry
    double CR5C_width[1]={0.253};                                                                     // the number of strips with equal pitches & width of the corresponding group of strips
    double CR5C_group[1]={1024};
    // For CR6 the numbers are final and should not be changed
    double CR6C_width[14]={0.38,0.32,0.27,0.23,0.17,0.18,0.22,0.25,0.29,0.33,0.37,0.41,0.46,0.51};    // the number of strips with equal pitches & width of the corresponding group of strips
    double CR6C_group[14]={32,32,32,32,704,64,32,32,32,32,32,32,32,32};

    int MxGrpSize = 14; // the max number of entries in CRC_group array
    CRCGROUP.resize(3); //3 regions
    CRCWIDTH.resize(3);

    for (int i = 0; i <3 ; ++i) {
    	CRCGROUP[i].resize(MxGrpSize);
    	CRCWIDTH[i].resize(MxGrpSize);
    }

    for(int j =0; j<13; j++)
    { // region index  0 is CR4
    	CRCGROUP[0][j] = CR4C_group[j];
    	CRCWIDTH[0][j] = CR4C_width[j];
    }
    for(int j =0; j<1; j++)
    { // region index  1 is CR5
        CRCGROUP[1][j] = CR5C_group[j];
        CRCWIDTH[1][j] = CR5C_width[j];
    }
    for(int j =0; j<14; j++)
    { // region index  2 is CR6
        CRCGROUP[2][j] = CR6C_group[j];
        CRCWIDTH[2][j] = CR6C_width[j];
    }

    CRCEDGE1.resize(3); //3 regions
	for (int i = 0; i <3 ; ++i)
		CRCEDGE1[i].resize(3);
	CRCEDGE2.resize(3); //3 regions
	for (int i = 0; i <3 ; ++i)
		CRCEDGE2[i].resize(3);

    double CEdge1[] = {toRadians(30.52),toRadians(270.52),toRadians(150.52)};
    double CEdge2[] = {toRadians(149.48),toRadians(29.48),toRadians(269.48)};

    // Assume the edge boundaries are the same for all regions until get final geometry
    for (int i = 0; i <3 ; ++i)
    {
    	for (int j = 0; j <3 ; ++j)
    	{
			CRCEDGE1[i][j] = CEdge1[j];
			CRCEDGE2[i][j] = CEdge2[j];
		}
    }

}

vector<double> bmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{

	double w_i = 25; //ionization potential assumed to be 25 eV
	int Nel = (int) (1e6*Edep/w_i);

	// the return vector is always in pairs the first index is the strip number, the second is the Edep on the strip
	vector<double> strip_id;

	int strip = -1;

	for(int iel=0;iel<Nel;iel++)
	{
		int num_region = (int) (layer+1)/2 - 1; 	// region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6;

		double sigma =0;
		if(layer%2==0)
		{// C layer
			double z_i = CRZZMIN[num_region]+CRZOFFSET[num_region]; 						 // fiducial z-profile lower limit
			double z_f = CRZZMIN[num_region]+CRZOFFSET[num_region] + CRZLENGTH[num_region];  // fiducial z-profile upper limit

			if(z>=z_i && z<=z_f)
			{
				sigma = getSigmaLongit(layer, x, y); //  longitudinal diffusion

				double z_d = (double) (G4RandGauss::shoot(z,sigma));

				strip = getCStrip(layer, z_d);

			}
		}

		if(layer%2==1)
		{// Z layer
			sigma = getSigmaAzimuth(layer, x, y); //   diffusion  taking into account the Lorentz angle

			double Delta_rad = sqrt(x*x+y*y)-CRZRADIUS[num_region]+hStrip2Det;

			double phi = atan2(y,x);
			double phi_d = phi + ((G4RandGauss::shoot(0,sigma))/cos(ThetaL)-Delta_rad*tan(ThetaL))/CRZRADIUS[num_region];

			strip = getZStrip(layer, phi_d);

		}

		if(strip != -1)
		{
			for(int istrip=0;istrip< (int) (strip_id.size()/2);istrip++)
			{
				if(strip_id[2*istrip]==strip)
				{// already hit strip - add Edep
					strip_id[2*istrip+1]=strip_id[2*istrip+1]+1./((double) Nel); // no gain fluctuation ;
					strip = -1; // not to use it anymore
				}
			}
			if(strip > -1) { // new strip
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
	if(Nel==0)
    { // Nel=0, consider the Edep is 0
        strip_id.push_back(-1);
        strip_id.push_back(1);
    }



    return strip_id;
}
double bmt_strip::toRadians(double angleDegrees) {
	return Pi*angleDegrees/180.;
}

/**
 *
 * param layer
 * param x x-coordinate of the hit in the lab frame
 * param y y-coordinate of the hit in the lab frame
 * return the sigma along the beam direction (longitudinal)
 */
double bmt_strip::getSigmaLongit(int layer, double x, double y)
{ // sigma for C-detector

	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	double sigma = SigmaDrift*sqrt((sqrt(x*x+y*y)-CRCRADIUS[num_region]+hStrip2Det)/hDrift);

	return sigma;
}

/**
 *
 * param layer
 * param x x-coordinate of the hit in the lab frame
 * param y y-coordinate of the hit in the lab frame
 * return the sigma along in the azimuth direction taking the Lorentz angle into account
 */
double bmt_strip::getSigmaAzimuth(int layer, double x, double y)
{ // sigma for Z-detectors

	int num_region = (int) (layer+1)/2 - 1; ///< region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	double sigma = SigmaDrift*sqrt((sqrt(x*x+y*y)-CRZRADIUS[num_region]+hStrip2Det)/hDrift/cos(ThetaL));

	return sigma;

}

/**
 *
 * param layer the layer 1...6
 * param angle the position angle of the hit in the Z detector
 * return the Z strip as a function of azimuthal angle
 */
int bmt_strip::getZStrip(int layer, double angle)
{
	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	int num_detector =isInSector( layer,  angle) - 1;
	if(num_detector==-1)
		return -1;

	if(angle<0)
		angle+=2*Pi; // from 0 to 2Pi

	if(num_detector==1) {
		double angle_f=CRCEDGE1[num_region][1]+(CRCXPOS[num_region]+CRCLENGTH[num_region])/CRCRADIUS[num_region] - 2*Pi;
		if(angle>=0 && angle<=angle_f)
			angle+=2*Pi;
	}
	double strip_calc = ( (angle-CRZEDGE1[num_region][num_detector])*CRZRADIUS[num_region]-CRZXPOS[num_region]-CRZWIDTH[num_region]/2.)/(CRZWIDTH[num_region]+CRZSPACING[num_region]);
	strip_calc = (int)(round(strip_calc) );
	int strip_num = (int) floor(strip_calc);

	int value = strip_num + 1;

	if(value<1 || value>CRZNSTRIPS[num_region])
		value = -1;

	return value;
}
/**
 *
 * param sector
 * param layer
 * param trk_z the track z position of intersection with the C-detector
 * return the C-strip
 */
int bmt_strip::getCStrip(int layer, double trk_z) {

	int strip_group = 0;
	int ClosestStp = -1;
	int StartStrip = 0;
	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6
	double Zb=0;
	double Z0 =0;
	double z =0;
	// get group
	int len = CRCGROUP[num_region].size();
	double Z_lowBound[len];
	double Z_uppBound[len];

	Z_lowBound[0] = CRCZMIN[num_region]+CRCOFFSET[num_region]+CRCWIDTH[num_region][0]/2.;
	Z_uppBound[0] = Z_lowBound[0]
					   + (CRCGROUP[num_region][0]-1)*(CRCWIDTH[num_region][0]/2.+ CRCSPACING[num_region])
					   + CRCWIDTH[num_region][0]/2;
	for(int i =1; i< len; i++)
	{
		Z_lowBound[i] = Z_uppBound[i-1] + CRCWIDTH[num_region][i-1]/2. ;
		Z_uppBound[i] = Z_lowBound[i] + (CRCGROUP[num_region][i]-1)*(CRCWIDTH[num_region][i]/2.+ CRCSPACING[num_region]);
		cout<<"  new bound "<<Z_lowBound[i]<<endl;
		if(trk_z>=Z_lowBound[i] && trk_z<=Z_uppBound[i]) {
			strip_group = i;
			cout<<"  new group "<<strip_group<<endl;
			break;
		}
	}

	for(int group =0; group< len; group++)
	{
		double zi= CRCZMIN[num_region]+CRCOFFSET[num_region];
		double z0 = CRCWIDTH[num_region][group]/2.;
		double zb = CRCGROUP[num_region][group]*(CRCWIDTH[num_region][group] + CRCSPACING[num_region]) ;

		z = trk_z - zi;
		z=round(z); // edge effect fix
		Z0+=z0;
		Zb=zb+Z0;

		if(z>=Z0 && z<Zb)
		{
			strip_group = group;
			if(strip_group>0) {
				for(int s =0; s<strip_group; s++)
					StartStrip+=CRCGROUP[num_region][strip_group];
			}
			break;
		}

	}

	int cClosestStrip = (int) (floor(((z-Z0)/(CRCWIDTH[num_region][strip_group] + CRCSPACING[num_region])))+0.5)+StartStrip;
	cout<<" strip_group "<<strip_group<<" c strip "<<cClosestStrip<<" start "<<StartStrip<<endl;
	for(int i =0; i< len; i++)
	{
		if(CRCGROUP[num_region][i]==0)
			break;
		int group =i;
		double zi= CRCZMIN[num_region]+CRCOFFSET[num_region];
		double z0 = CRCWIDTH[num_region][group]/2.;
		double zb = CRCGROUP[num_region][group]*(CRCWIDTH[num_region][group] + CRCSPACING[num_region]) ;

		double z = trk_z - zi;
		z=round(z); // edge effect fix
		Z0+=z0;
		Zb=zb+Z0;

		if(z>=Z0 && z<Zb)
		{
			strip_group = group;
			int min_strip = 1;
			for(int g =0; g<strip_group; g++)
				min_strip+=CRCGROUP[num_region][g];
			int max_strip = min_strip + CRCGROUP[num_region][strip_group];

			double StripDiffMin = CRCLENGTH[num_region];

			for(int nCstrpNb = min_strip; nCstrpNb<=max_strip; nCstrpNb++)
			{
				double zstp = CRCStrip_GetZ(layer, nCstrpNb ) -zi; //  c strip
				zstp = ceil(zstp*100000)/100000; // rounding fix

				double StripDiffCalc = abs(z-zstp);
				if(StripDiffCalc<StripDiffMin)
				{
					StripDiffMin = StripDiffCalc;
					ClosestStp = nCstrpNb;

				}
				else
				{
					break;
				}
			}
		}
		Z0+=zb;
	}


	if(ClosestStp<1 || ClosestStp>CRCNSTRIPS[num_region])
		ClosestStp = -1;
	cout<<" Closeststp "<<ClosestStp<<endl;
	return ClosestStp;
}

/**
 * param layer the hit layer
 * param strip the hit strip
 * return the z position in mm for the C-detectors
 */
double bmt_strip::CRCStrip_GetZ(int layer, int strip)
{
	int num_strip = strip - 1;     			// index of the strip (starts at 0)
	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	//For CR6C, this function returns the Z position of the strip center
	int group=0;
	int limit=CRCGROUP[num_region][group];
	double zc=CRCZMIN[num_region]+CRCOFFSET[num_region]+CRCWIDTH[num_region][group]/2.;

	if (num_strip>0){
	  for (int j=1;j<num_strip+1;j++)
	  {
		zc+=CRCWIDTH[num_region][group]/2.;
		if (j>=limit)
		{ //test if we change the width
			group++;
			limit+=CRCGROUP[num_region][group];
		}
		zc+=CRCWIDTH[num_region][group]/2.+CRCSPACING[num_region];
	  }
	}

	return zc; //in mm
}
/**
*
* param sector the sector in CLAS12 1...3
* param layer the layer 1...6
* param strip the strip number (starts at 1)
* return the angle to localize the  center of strip
*/
double bmt_strip::CRZStrip_GetPhi(int sector, int layer, int strip){
	// Sector = num_detector + 1;
	// num_detector = 0 (region A), 1 (region B), 2, (region C)
	//For CRZ, this function returns the angle to localize the  center of strip "num_strip" for the "num_detector"
	int num_detector = getDetectorIndex(sector); 			// index of the detector (0...2): sector 1 corresponds to detector B, 2 to A,  3 to C in the geometry created by GEMC
	int num_strip = strip - 1;     							// index of the strip (starts at 0)
	int num_region = (int) (layer+1)/2 - 1; 				// region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	double angle=CRZEDGE1[num_region][num_detector]+(CRZXPOS[num_region]+(CRZWIDTH[num_region]/2.+num_strip*(CRZWIDTH[num_region]+CRZSPACING[num_region])))/CRZRADIUS[num_region];
	if (angle>2*Pi) angle-=2*Pi;
	return angle; //in rad
}

double bmt_strip::getEnergyFraction(double z0, double z, double sigma){
	double pdf_gaussian = (1./(sigma*sqrt(2*Pi)))* exp( -0.5*((z-z0)/sigma)*((z-z0)/sigma) );
	return pdf_gaussian;
}

int bmt_strip::getDetectorIndex(int sector) {
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

int bmt_strip::isInSector(int layer, double angle) {

	int num_region = (int) (layer+1)/2 - 1; // region index (0...2) 0=layers 1&2, 1=layers 3&4, 2=layers 5&6

	if(angle<0)
		angle+=2*Pi; // from 0 to 2Pi
	double angle_pr = angle + 2*Pi;

	double angle_i = 0; // first angular boundary init
	double angle_f = 0; // second angular boundary for detector A, B, or C init
	int num_detector = -1;

	for(int i = 0; i<3; i++) {

		angle_i=CRCEDGE1[num_region][i]+CRCXPOS[num_region]/CRCRADIUS[num_region];
		angle_f=CRCEDGE1[num_region][i]+(CRCXPOS[num_region]+CRCLENGTH[num_region])/CRCRADIUS[num_region];

		if( (angle>=angle_i && angle<=angle_f) || (angle_pr>=angle_i && angle_pr<=angle_f) )
			num_detector=i;
	}
	return num_detector + 1;
}
