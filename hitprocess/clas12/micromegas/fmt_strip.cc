// gemc headers
#include "fmt_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>

void fmt_strip::fill_infos()
{
  // all dimensions are in mm
  
  Pi  = 3.14159265358;
  //  interlayer      = (0.5+0.1+0.015+0.128+0.030+2.500)*2.; // distance between 2 layers of a superlayer // not usefull anymore, 07/27/2015 (Frederic Georges)
  //  intersuperlayer = 20.0;      // distance between 2 superlayers
  //  intersuperlayer = 21.0;      // distance between 2 superlayers // modified on 5/20/2015 to match new geometry (R.De Vita)
  intersuperlayer = 10.5;      // distance between 2 superlayers // modified on 7/27/2015 to match new geometry (Frederic Georges)
  //  pitch           = 0.500;     // pitch of the strips
  pitch           = 0.525;     // pitch of the strips // modified on 7/27/2015 to match new geometry (Frederic Georges)
  hDrift          = 5.0;
  //hStrip2Det      = hDrift/2.;
  sigma_td_max    = 0.01; // very small transverse diffusion because of B field
  w_i             = 25.0;
  
  //R_min = 12.900;
  //R_max = 215.0;  // outer radius of strip part
  R_min = 45.500;   // modified on 7/27/2015 to match new geometry (Frederic Georges)
  R_max = 185.400;  // outer radius of strip part // modified on 7/27/2015 to match new geometry (Frederic Georges)
  //  Z_1stlayer = 295.-0.5-0.1-0.015-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap)
  //Z_1stlayer = 305.25-0.5-0.1-0.015-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap) // modified on 5/20/2015 to match new geometry (R.De Vita)
  Z_1stlayer = 305.250-0.005-0.200-2.000-0.200-0.005-0.050-0.020-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap) // modified on 7/27/2015 to match new geometry (Frederic Georges)
  
  // z of the upstream part of the layer
  /*
  Z0.push_back(Z_1stlayer);
  Z0.push_back(Z0[0]+interlayer);
  Z0.push_back(Z_1stlayer+intersuperlayer);
  Z0.push_back(Z0[2]+interlayer);
  Z0.push_back(Z_1stlayer+2.*intersuperlayer);
  Z0.push_back(Z0[4]+interlayer);
  */
  // modified on 7/27/2015 to match new geometry (Frederic Georges)
  Z0.push_back(Z_1stlayer);
  Z0.push_back(Z0[0]+intersuperlayer);
  Z0.push_back(Z0[1]+intersuperlayer);
  Z0.push_back(Z0[2]+intersuperlayer);
  Z0.push_back(Z0[3]+intersuperlayer);
  Z0.push_back(Z0[4]+intersuperlayer);
  
  // angles of each layer
  /*
  alpha.push_back(0);
  alpha.push_back(Pi/2.);
  alpha.push_back(Pi/3.);
  alpha.push_back(Pi/2+Pi/3.);
  alpha.push_back(2.*Pi/3.);
  alpha.push_back(2.*Pi/3+Pi/2.);
  */
  // modified on 7/27/2015 to match new geometry (Frederic Georges)
  alpha.push_back(19.*Pi/180.);
  alpha.push_back(alpha[0]+Pi/3.);
  alpha.push_back(alpha[0]+2.*Pi/3.);
  alpha.push_back(alpha[0]+Pi);
  alpha.push_back(alpha[0]+4.*Pi/3.);
  alpha.push_back(alpha[0]+5.*Pi/3.);
  
  // Number of strips and pixels
  N_str = 1024; //16 connectors * 64 strips = 1024 strips for each fmt
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
	  x_real = (double) (G4RandGauss::shoot(x*cos(alpha[layer])+y*sin(alpha[layer]),sigma_td));
	  y_real = (double) (G4RandGauss::shoot(y*cos(alpha[layer])-x*sin(alpha[layer]),sigma_td));
	  
	  if(y_real > -84.1 && y_real < 84.1 && x_real < 0){ // R_max - 3*64*0.525 - 0.5 = 84.1;
		ClosestStrip = (int) (floor((84.1-y_real)/pitch)+1);
	  	}
	  else if(y_real < -84.1 && y_real > -184.9){ // R_max - 3*64*0.525 - 0.5 = 84.1; R_max - 0.5 = 184.9
	  	ClosestStrip = (int) (floor((-y_real-84.1)/pitch)+1) + 320; // 5*64 = 320
		}
	  else if(y_real > -84.1 && y_real < 84.1 && x_real > 0){ // R_max - 3*64*0.525 - 0.5 = 84.1;
		ClosestStrip = (int) (floor((y_real+84.1)/pitch)+1) + 512;	 // (5+3)*64 = 512
		}
	  else if(y_real > 84.1 && y_real < 184.9){ // R_max - 3*64*0.525 - 0.5 = 84.1; R_max - 0.5 = 184.9
		  ClosestStrip = (int) (floor((y_real-84.1)/pitch)+1) + 832; // (5+3+5)*64 = 832
		}
	
	  if(sqrt(x*x+y*y)<R_max && sqrt(x*x+y*y)>R_min && ClosestStrip>=1 && ClosestStrip<=N_str)
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











