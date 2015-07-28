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
 interlayer 	= 16.0;
 pitchZ      	= 0.540;
 pitchC      	= 0.270;
 Nsector    	= 3;
 hDrift         = 3.0;
 hStrip2Det 	= hDrift/2.;
 sigma_td_max   = 0.4;
 theta_L        = Pi*20./180.;
 w_i            = 25.0;

 DZ_inLength = 7.5;
 DZ_inWidth  = 7.5;

 // z of the upstream part of the layer
 Z0.push_back(-127.000);  Z0.push_back(-127.000);
 Z0.push_back(-148.000); Z0.push_back(-148.000);
 Z0.push_back(-172.500); Z0.push_back(-172.500);
 
// total z length of the layer of the layer
 DZ.push_back(370.5);  DZ.push_back(370.5);
 DZ.push_back(420.1); DZ.push_back(420.1);
 DZ.push_back(444.6); DZ.push_back(444.6);

 // radii of layers
 R.push_back(146.900); R.push_back(R[0]+interlayer);
 R.push_back(178.900); R.push_back(R[2]+interlayer);
 R.push_back(210.900); R.push_back(R[4]+interlayer);

 // Number of strips (depends on radius, dead zones, and pitch!)
 for(int i=0;i<6;i++)
   { // layer 0,2,4 are Z detectors, i.e. measure phi (strips along z)
     if((i%2)==0) Nstrips.push_back((int) ((2.*Pi*R[i]/((double) Nsector)-2.*DZ_inLength)/pitchZ));
     if((i%2)==1) Nstrips.push_back((int) ((DZ[i]-2.*DZ_inWidth)/pitchC));
   }
 
 // mid angle of the sector
 MidTile.push_back(0);     MidTile.push_back(0);
 MidTile.push_back(0);     MidTile.push_back(0);
 MidTile.push_back(0);     MidTile.push_back(0);
}

vector<double>  bmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
  // the return vector is always in pairs. 
  // The first number is the ID, 
  // the second number is the sharing percentage
  vector<double> strip_id;
  
 // number of electrons (Nt)
 Nel = (int) (1e6*Edep/w_i);
 
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
 int ClosestStrip=0;
 
 // now compute the sigma of the (transverse) dispersion for this interaction
 if((layer%2)==1) sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/hDrift); // "C" det, transverse diffusion grows with square root of distance
 else sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/(cos(theta_L)*hDrift)); // same, but "Z" detectors, so Lorentz angle makes drift distance longer by 1./cos(theta_L) . Means sigma_td can be larger than sigma_td_max

 if(Nel>0)
   {
     for(int iel=0;iel<Nel;iel++)
       { // loop over (total) electrons
	 if((layer%2)==1) 
	   { // this for "C" layers, i.e. measuring z
	     z_real = (double) (G4RandGauss::shoot(z,sigma_td));
	     ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth)/pitchC))+0.5);
	   }
	 if((layer%2)==0)
	   { // this for "Z" layers, i.e. measuring phi
	     phi_real = phi + ((G4RandGauss::shoot(0,sigma_td))/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer]; // the sign of the 2nd term (Lorentz angle) should be a "-" as the B field is along +z and the MM are convex (and electrons are negatively charged)
	     ClosestStrip = (int) (floor(((R[layer]/pitchZ)*(phi_real-phiij+Pi/Nsector-DZ_inLength/R[layer]))+0.5));
	   }
	 if(ClosestStrip>=0 && ClosestStrip<=Nstrips[layer] && z>=Z0[layer]+DZ_inWidth && z<=Z0[layer]+DZ[layer]-DZ_inWidth && phi>=phiij-Pi/Nsector+DZ_inLength/R[layer] && phi<=phiij+Pi/Nsector-DZ_inLength/R[layer])
	   { // strip is in the acceptance, check if new strip or not
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
	   { // not in the acceptance
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

