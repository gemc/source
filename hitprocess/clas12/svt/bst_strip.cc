// gemc headers
#include "bst_strip.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

void bst_strip::fill_infos()
{
	// ccdb path : "/geometry/cvt/svt/";
	
	intersensors = (0.11 + 0.002)*mm;  // gap between sensors + microgap
	// intersensors = /geometry/cvt/svt/svt/microGapLen
	alpha        = 3.0*deg;            // max angle of the strips
	//alpha = /geometry/cvt/svt/svt/stereoAngle
	pitch        = 0.156*mm;           // pitch of the strips
	//pitch = /geometry/cvt/svt/svt/readoutPitch

	//========================//
	// other variables needed //
	//========================//
	StripStart 	= 0.048*mm;					 // the first sensitive strip starts at stripStart + 1/2 of the pitch
	//stripStart = /geometry/cvt/svt/svt/stripStart

	// number of cards by sector for each layer
	NSensors.push_back(3);   NSensors.push_back(3);
	NSensors.push_back(3);   NSensors.push_back(3);
	NSensors.push_back(3);   NSensors.push_back(3);
	NSensors.push_back(3);   NSensors.push_back(3);
	
	DZ_inWidth 		= 0.984*mm;  	 // size of dead zone in the direction of the length of the module
	//DZ_inWidth = /geometry/cvt/svt/svt/deadZnWid
	DZ_inLength   	= 0.835*mm;  	 // size of dead zone in the direction of the width of the module
	//DZ_inLength = /geometry/cvt/svt/svt/deadZnLen
	SensorLength 	= 111.625*mm;	 // length of 1 Sensor
	//SensorLength = /geometry/cvt/svt/svt/physSenLen
	SensorWidth  	= 42.0*mm; 		 // width 1 Sensor
	//SensorWidth = /geometry/cvt/svt/svt/physSenWid
	
	// Number of strips -1
	Nstrips  = (int) floor((SensorWidth-2.0*DZ_inWidth)/pitch) - 1;
	
}


vector<double> bst_strip::FindStrip(int layer, int sector, int isens, G4ThreeVector Lxyz)
{
	
	// the return vector is always in pairs.
	// The first number is the ID,
	// the second number is the sharing percentage
	// layer and sector are the indexes (i.e. layer is the actual layern-1)
	vector<double> strip_id;
	
	int StripHit = -1;
	double minDist = 999;  // distance between actual x and x of the k-strip at the actual z
	double dalpha = alpha/((float)Nstrips);
	double SensitiveSensorWidth = (SensorWidth - 2.0*DZ_inWidth); // sensitive width
	
	// local x position in the sensor: need to add 1/2 active area
	//double lx = Lxyz.x() + SensorWidth/2.0 - DZ_inWidth;
	double lx = Lxyz.x() + SensitiveSensorWidth/2.0;
	// local z position in the module: need to add 1/2 active area + 1(2) full sensor lengths if on sensor 2(3)
	double lz = Lxyz.z() + SensorLength/2.0 - DZ_inLength + (isens-1)*(SensorLength + intersensors) ;
	
	// particle must be in the z active area of the sensor
	if( fabs(Lxyz.z()) < (SensorLength/2.0 - DZ_inLength) && fabs(Lxyz.x()) < SensitiveSensorWidth/2.0 )
	{
		
		// looping over all strips to find the closest
		// to the point
		for(int k=0; k<Nstrips; k++)
		{
			// angle of the k-strip
			double alpha_k = k*dalpha;
			
			// strip equation is z = mz + intcp
			double intcp = 0;
			double     m = 0;
			
			if(layer%2 == 0)
			{
				// intercept: for layer A it starts from the top of the active area of the sensor in sector 1
				// which is at positive local x
				// the first strip starts at StripStart+1/2Pitch; hence in this frame the first strip intercept is SensitiveSensorWidth-(StripStart+1/2Pitch)
				intcp = SensitiveSensorWidth - StripStart - (k+0.5)*pitch ;
				m     = -tan(alpha_k/rad);
			}
			
			if(layer%2 == 1)
			{
				// intercept: for layer B it starts from the bottom of the active area of the sensor in sector 1
				// which is at negative local x
				// the first strip starts at StripStart+1/2Pitch
				intcp = StripStart + (k+0.5)*pitch ;
				m     = tan(alpha_k/rad);
			}
			
			// x position of the strip according to the active local z position
			double stripx = intcp + m*lz;
			
			if(fabs(lx - stripx) < fabs(minDist))
			{
				minDist = lx - stripx;
				StripHit = k+1;
			}
			
		}
	}
	
	
	
	if(StripHit != -1)
	{
		// correcting for z positioning inside the module
		double dpitch = pitch + lz*tan(dalpha/rad);

		// for even layers x is proportional to strip number.
		// for odd layers x is inversly proportional to strip number
		int moduleDirection = 1;
		if(layer%2 == 0)
			moduleDirection = -1;

		// one strip only, all energy to it
		if(fabs(minDist)<=dpitch/4.0)
		{
			strip_id.push_back(StripHit);
			strip_id.push_back(1);
		}
		// two hits, on the right of the strip (left if odd layers)
		// 10% loss due to capacitance between strip and backplane
		if(minDist>dpitch/4.0)
		{
			strip_id.push_back(StripHit);
			strip_id.push_back(0.45);
			if(StripHit + moduleDirection != 0)
			{
				strip_id.push_back(StripHit + moduleDirection);
				strip_id.push_back(0.45);
			}
			//if(StripHit == 1) cout << "charged shared with " << StripHit + moduleDirection << endl;
		}
		// two hits, on the left of the strip (right if odd layers)
		// 10% loss due to capacitance between strip and backplane
		if(minDist<-dpitch/4.0)
		{
			strip_id.push_back(StripHit);
			strip_id.push_back(0.45);
			if(StripHit - moduleDirection != 0)
			{
				strip_id.push_back(StripHit - moduleDirection);
				strip_id.push_back(0.45);
			}
			//if(StripHit == 1) cout << "charged shared with " << StripHit + moduleDirection << endl;
		}
	}
	else
	{
		strip_id.push_back(-5);
		strip_id.push_back(1);
	}
	return strip_id;
	
}


