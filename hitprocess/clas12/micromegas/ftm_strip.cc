// gemc headers
#include "ftm_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>

void ftm_strip::fill_infos()
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
  
  // z of the upstream part of the layer
  Z0.push_back(Z_1stlayer);
  Z0.push_back(Z0[0]+interlayer);
  Z0.push_back(Z_1stlayer+intersuperlayer);
  Z0.push_back(Z0[2]+interlayer);
  
  // Number of strips
  Nstrips = (int) floor(2.*(Rmax+Rmin)/pitch);
}

vector<double> ftm_strip::FindStrip(int layer, double x, double y, double z, double Edep)
{
  // the return vector is always in pairs. 
  // The first number is the ID, 
  // the second number is the sharing percentage
  vector<double> strip_id;
  // number of electrons (Nt)
  Nel = (int) (1e6*Edep/w_i);
  if(fabs(z-Z0[layer])>(hDrift+0.2)) cout << "Warning! z position of the FTM hit is not in the sensitive volume: " << z-Z0[layer]<< endl;
  sigma_td = sigma_td_max* sqrt(fabs(z-Z0[layer])/hDrift); // expression without Lorentz angle
  
  int ClosestStrip=0;
  if(Nel>0)
    {
      for(int iel=0;iel<Nel;iel++)
	{ // loop over (total) electrons
	  if(layer%2==0)
	    {
	      x_real = (double) (G4RandGauss::shoot(x,sigma_td));
	      y_real = y;
	      if(-Rmin<x_real && x_real<Rmin && y_real<0) ClosestStrip = (int) (floor(2.0*Rmax/pitch+(x_real+Rmin)/pitch+0.5));
	      else ClosestStrip = (int) (floor((x_real+Rmax)/pitch+0.5));
	    }
	  if(layer%2==1)
	    {
	      x_real = x;
	      y_real = (double) (G4RandGauss::shoot(y,sigma_td));
	      if(-Rmin<y_real && y_real<Rmin && x_real<0) ClosestStrip = (int) (floor(2.0*Rmax/pitch+(y_real+Rmin)/pitch+0.5));
	      else ClosestStrip = (int) (floor((y_real+Rmax)/pitch+0.5));
	    }
	  
	  if(sqrt(x_real*x_real+y_real*y_real)<Rmax && sqrt(x_real*x_real+y_real*y_real)>Rmin && ClosestStrip>=0 && ClosestStrip<=Nstrips)
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

