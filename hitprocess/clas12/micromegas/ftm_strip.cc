// gemc headers
#include "ftm_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>

/*
void ftm_strip::fill_infos(detector Detector)
{
  // all dimensions are in mm
 
  Pi  = 3.14159265358;
  interlayer      = (0.3/2.+0.1+0.015)*2.; // distance between 2 layers of a superlayer
  intersuperlayer = 20.0;      // distance between 2 superlayers
  
  pitch           = 0.500;                         // pitch of the strips
  hDrift          = 5.35;
  hStrip2Det      = hDrift/2.;
  sigma_td_max    = 0.0; // very small transverse diffusion (temporary)
  w_i             = 25.0;
  
  Rmin = 65.0;                 // inner radius of disks
  Rmax = 142.0;  // outer radius of disks
  Z_1stlayer = 1773.-0.3/2.-0.1-0.015; // z position of the 1st layer : epoxy center - half epoxy Dz - PCB Dz - strips Dz
  
    cout << Detector.name << " " << Detector.dimensions.size();
    for(int i=0; i<Detector.dimensions.size(); i++) cout << " " << Detector.dimensions[i];
    cout << endl;
  // z of the upstream part of the layer
  Z0.push_back(Z_1stlayer);
  Z0.push_back(Z0[0]+interlayer);
  Z0.push_back(Z_1stlayer+intersuperlayer);
  Z0.push_back(Z0[2]+interlayer);
  
  // Number of strips
  Nstrips = (int) floor(2.*(Rmax+Rmin)/pitch);
}
*/
    
vector<double> ftm_strip::FindStrip(int layer, double x, double y, double z, double Edep, detector Detector,ftmConstants ftmcc)
{
    // the return vector is always in pairs.
    // The first number is the ID,
    // the second number is the sharing percentage
    vector<double> strip_id;
    // number of electrons (Nt)
    int Nel = (int) (1e6*Edep/ftmcc.w_i);
    
    // get layer position and drift distance from the detector dimensions
    double z0 = -Detector.dimensions[2];
    double hdrift = Detector.dimensions[2]*2;
//    cout << Detector.name << " " << Detector.dimensions.size();
//    for(int i=0; i<Detector.dimensions.size(); i++) cout << " " << Detector.dimensions[i];
//    cout << endl;
//    cout << " z0 " <<z0 << " drift " << hdrift << " z " << z << endl;

    // old if(fabs(z-Z0[layer])>(hDrift+0.2)) cout << "Warning! z position of the FTM hit is not in the sensitive volume: " << z-Z0[layer]<< endl;
    if(fabs(z-z0)>(hdrift+0.2)) cout << "Warning! z position of the FTM hit is not in the sensitive volume: " << z << " " << z0<< endl;
    
    double sigma_td = ftmcc.sigma_td_max* sqrt(fabs(z-z0)/hdrift); // expression without Lorentz angle
//    cout << Nel << endl;
    int ClosestStrip=0;
    if(Nel>0)
    {
      for(int iel=0;iel<Nel;iel++)
	{ // loop over (total) electrons
        double x_real=x;
        double y_real=y;
	  if(layer%2==0)
	    {
          x_real = (double) (G4RandGauss::shoot(x,sigma_td));
          y_real = y;
	      if(-ftmcc.rmin<x_real && x_real<ftmcc.rmin && y_real<0) ClosestStrip = (int) (floor(2.0*ftmcc.rmax/ftmcc.pitch+(x_real+ftmcc.rmin)/ftmcc.pitch+0.5));
	      else ClosestStrip = (int) (floor((x_real+ftmcc.rmax)/ftmcc.pitch+0.5));
	    }
	  if(layer%2==1)
	    {
	      x_real = x;
	      y_real = (double) (G4RandGauss::shoot(y,sigma_td));
	      if(-ftmcc.rmin<y_real && y_real<ftmcc.rmin && x_real<0) ClosestStrip = (int) (floor(2.0*ftmcc.rmax/ftmcc.pitch+(y_real+ftmcc.rmin)/ftmcc.pitch+0.5));
	      else ClosestStrip = (int) (floor((y_real+ftmcc.rmax)/ftmcc.pitch+0.5));
	    }
	  
	  if(sqrt(x_real*x_real+y_real*y_real)<ftmcc.rmax && sqrt(x_real*x_real+y_real*y_real)>ftmcc.rmin && ClosestStrip>=0 && ClosestStrip<=ftmcc.nstrips)
	    { // strip is in the acceptance
	      for(int istrip=0;istrip< (int) (strip_id.size()/2);istrip++)
		{
		  if(strip_id[2*istrip]==ClosestStrip)
		    {// already hit strip - add Edep
		      strip_id[2*istrip+1]=strip_id[2*istrip+1]+1./((double) Nel); // no gain fluctuation yet
		      ClosestStrip=-1; // not to use it anymore
		    }
		}
	      if(ClosestStrip>-1)
		{ // this is a new strip
		  strip_id.push_back(ClosestStrip);
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

