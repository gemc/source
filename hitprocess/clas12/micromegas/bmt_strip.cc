// gemc headers
#include "bmt_strip.h"
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
    ThetaL = toRadians(20.);

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
    for (int i = 0; i <3 ; ++i)	{
    	for (int j = 0; j <3 ; ++j){
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
    double width[1]={1024};
    // For CR6 the numbers are final and should not be changed
    double CR6C_width[14]={0.38,0.32,0.27,0.23,0.17,0.18,0.22,0.25,0.29,0.33,0.37,0.41,0.46,0.51};    // the number of strips with equal pitches & width of the corresponding group of strips
    double CR6C_group[14]={32,32,32,32,704,64,32,32,32,32,32,32,32,32};

    int MxGrpSize = 14; // the max number of entries in CRC_group array
    CRCGROUP.resize(3); //3 regions
    CRCWIDTH.resize(3);

    for (int i = 0; i <3 ; ++i)
    	CRCGROUP[i].resize(MxGrpSize);
    for(int j =0; j<sizeof(CR4C_group) / sizeof(CR4C_group[0]); j++) { // region index  0 is CR4
    	CRCGROUP[0][j] = CR4C_group[j];
    	CRCWIDTH[0][j] = CR4C_width[j];
    }
    for(int j =0; j<sizeof(CR5C_group) / sizeof(CR5C_group[0]); j++) { // region index  1 is CR5
        CRCGROUP[1][j] = CR5C_group[j];
        CRCWIDTH[1][j] = CR5C_width[j];
    }
    for(int j =0; j<sizeof(CR6C_group) / sizeof(CR6C_group[0]); j++) { // region index  2 is CR6
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
    for (int i = 0; i <3 ; ++i)	{
    	for (int j = 0; j <3 ; ++j){
			CRCEDGE1[i][j] = CEdge1[j];
			CRCEDGE2[i][j] = CEdge2[j];
		}
    }

}
/**
 * Method to get digi hits based on x,y,z position and E
 */
vector<double>  bmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
        // the return vector is always in pairs.
	vector<double> strip_id;
		{ // Nel=0, consider the Edep is 0
				strip_id.push_back(-1);
				strip_id.push_back(0);
		}

        return strip_id;
}


double bmt_strip::toRadians(double angleDegrees) {
	return Pi*angleDegrees/180.;
}
