// gemc headers
#include "fmt_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>

void fmt_strip::fill_infos()
{
  // all dimensions are in mm
  
  Pi  = 3.14159265358;
  interlayer      = (0.5+0.1+0.015+0.128+0.030+2.500)*2.; // distance between 2 layers of a superlayer
  //  intersuperlayer = 20.0;      // distance between 2 superlayers
  intersuperlayer = 21.0;      // distance between 2 superlayers // modified on 5/20/2015 to match new geometry (R.De Vita)
  pitch           = 0.500;     // pitch of the strips
  hDrift          = 5.0;
  hStrip2Det      = hDrift/2.;
  sigma_td_max    = 0.01; // very small transverse diffusion because of B field
  w_i             = 25.0;
  
  R_min = 12.900;
  R_max = 215.0;  // outer radius of strip part
  //  Z_1stlayer = 295.-0.5-0.1-0.015-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap)
  Z_1stlayer = 305.25-0.5-0.1-0.015-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap) // modified on 5/20/2015 to match new geometry (R.De Vita)
  
  // z of the upstream part of the layer
  Z0.push_back(Z_1stlayer);
  Z0.push_back(Z0[0]+interlayer);
  Z0.push_back(Z_1stlayer+intersuperlayer);
  Z0.push_back(Z0[2]+interlayer);
  Z0.push_back(Z_1stlayer+2.*intersuperlayer);
  Z0.push_back(Z0[4]+interlayer);
  
  // angles of each layer
  alpha.push_back(0);
  alpha.push_back(Pi/2.);
  alpha.push_back(Pi/3.);
  alpha.push_back(Pi/2+Pi/3.);
  alpha.push_back(2.*Pi/3.);
  alpha.push_back(2.*Pi/3+Pi/2.);
  
  // Number of strips and pixels
  N_str = (int) floor(2.*R_max/pitch);
}

vector<double> fmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
  // the return vector is always in pairs. 
  // The first number is the ID, 
  // the second number is the sharing percentage
  vector<double> strip_id;
  // number of electrons (Nt)
  Nel = (int) (1e6*Edep/w_i);
  if(z-Z0[layer]>0.5001*hDrift || z-Z0[layer]<-0.5001*hDrift) cout << "Warning! z position of the FMT hit is not in the sensitive volume: " << z-Z0[layer]<< endl;
  sigma_td = sigma_td_max* sqrt((z-Z0[layer])/hDrift); // expression without Lorentz angle
  
  int ClosestStrip=0;
  if(Nel>0)
    {
      for(int iel=0;iel<Nel;iel++)
	{ // loop over (total) electrons
	  u_real = (double) (G4RandGauss::shoot(y*cos(alpha[layer])-x*sin(alpha[layer]),sigma_td));
	  ClosestStrip = (int) (floor((u_real+R_max)/pitch));
	  
	  if(sqrt(x*x+y*y)<R_max && sqrt(x*x+y*y)>R_min && ClosestStrip>=0 && ClosestStrip<=N_str)
	    { // strip is in the acceptance (~)
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











