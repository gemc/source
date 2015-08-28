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
 interlayer 	= 15.0;
 Nsector    	= 3;
 hDrift         = 3.0;
 hStrip2Det 	= hDrift/2.;
 sigma_td_max   = 0.4;
 theta_L        = Pi*20./180.;
 w_i            = 25.0;


 interStripZ = 0.201;
 interStripC = 0.160;

 // pitch CR4Z
 pitchZ4 = 0.487; // and interpitch=0.201
 
 // pitch CR5Z
 pitchZ5 = 0.49; // none existing value. picked a rabdom value compatible with the geometry
 
 // pitch CR6Z
 pitchZ6 = 0.526; // interpitch=0.201 ; 0.526 = 0.201 + 0.325

 // pitch CR4C
 pitchC4.push_back(0.505); // 32 strips
 pitchC4.push_back(0.440); // 32 strips
 pitchC4.push_back(0.385); // 32 strips
 pitchC4.push_back(0.335); // 32 strips
 pitchC4.push_back(0.330); // 8 * 64 strips
 pitchC4.push_back(0.370); // 32 strips
 pitchC4.push_back(0.420); // 32 strips
 pitchC4.push_back(0.470); // 32 strips
 pitchC4.push_back(0.530); // 32 strips
 pitchC4.push_back(0.600); // 32 strips
 pitchC4.push_back(0.675); // 32 strips
 pitchC4.push_back(0.765); // 32 strips
 pitchC4.push_back(0.860); // 32 strips
 						   // 14 * 64
 
 // pitch CR5C
 // none existing value. picked a random value compatible with the geometry
 pitchC5.push_back(0.413); // 16 groupes de 64
 						   // 16 * 64
 
 // pitch CR6C
 pitchC6.push_back(0.54); // 32 strips
 pitchC6.push_back(0.48); // 32 strips
 pitchC6.push_back(0.43); // 32 strips
 pitchC6.push_back(0.39); // 32 strips
 pitchC6.push_back(0.33); // 11 * 64 strips
 pitchC6.push_back(0.34); // 32 strips
 pitchC6.push_back(0.38); // 32 strips
 pitchC6.push_back(0.41); // 32 strips
 pitchC6.push_back(0.45); // 32 strips
 pitchC6.push_back(0.49); // 32 strips
 pitchC6.push_back(0.53); // 32 strips
 pitchC6.push_back(0.57); // 32 strips
 pitchC6.push_back(0.62); // 32 strips
 pitchC6.push_back(0.67); // 32 strips
 						  // 18 * 64

 // z of the upstream part of the layer
 Z0.push_back(-127.820);  Z0.push_back(-127.820);
 Z0.push_back(-148.740); Z0.push_back(-148.740);
 Z0.push_back(-169.710); Z0.push_back(-169.710);
 
// total z length of the layer of the layer
 DZ.push_back(372.75);  DZ.push_back(372.75);
 DZ.push_back(423.99); DZ.push_back(423.99);
 DZ.push_back(444.96); DZ.push_back(444.96);

 // radii of layers
 R.push_back(145.731); R.push_back(R[0]+interlayer);
 R.push_back(175.731); R.push_back(R[2]+interlayer);
 R.push_back(205.731); R.push_back(R[4]+interlayer);

 // Number of strips (depends on radius, dead zones, and pitch!)
 for(int i=0;i<6;i++)
   { // layer 0,2,4 are Z detectors, i.e. measure phi (strips along z)
     if(i==0) Nstrips.push_back(14*64); // CR4C
     if(i==1) Nstrips.push_back(10*64); // CR4Z
	 if(i==2) Nstrips.push_back(11*64); // CR5Z
	 if(i==3) Nstrips.push_back(16*64); // CR5C
	 if(i==4) Nstrips.push_back(12*64); // CR6Z
	 if(i==5) Nstrips.push_back(18*64); // CR6C
   }
 
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

vector<double>  bmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
  // the return vector is always in pairs. 
  // The first number is the ID, 
  // the second number is the sharing percentage
  vector<double> strip_id;
  
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
 if(layer==0 || layer==3 || layer==5) sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/hDrift); // "C" det, transverse diffusion grows with square root of distance
 else sigma_td = sigma_td_max* sqrt((sqrt(x*x+y*y)-R[layer]+hStrip2Det)/(cos(theta_L)*hDrift)); // same, but "Z" detectors, so Lorentz angle makes drift distance longer by 1./cos(theta_L) . Means sigma_td can be larger than sigma_td_max


 if(Nel>0)
   {
     for(int iel=0;iel<Nel;iel++)
       { // loop over (total) electrons
	 
	 if(layer==0) 
	   { // this for "C4" layers, i.e. measuring z
	  	 z_real = (double) (G4RandGauss::shoot(z,sigma_td));
		 if(z_real-Z0[layer]-DZ_inWidth>0 && z_real-Z0[layer]-DZ_inWidth<=16.16)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth)/pitchC4[0]))+0.5);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>16.16 && z_real-Z0[layer]-DZ_inWidth<=30.24)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-16.16)/pitchC4[1]))+0.5+32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>30.24 && z_real-Z0[layer]-DZ_inWidth<=42.56)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-30.24)/pitchC4[2]))+0.5+2*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>42.56 && z_real-Z0[layer]-DZ_inWidth<=53.28)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-42.56)/pitchC4[3]))+0.5+3*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>53.28 && z_real-Z0[layer]-DZ_inWidth<=222.24)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-53.28)/pitchC4[4]))+0.5+4*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>222.24 && z_real-Z0[layer]-DZ_inWidth<=234.08)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-222.24)/pitchC4[5]))+0.5+20*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>234.08 && z_real-Z0[layer]-DZ_inWidth<=247.52)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-234.08)/pitchC4[6]))+0.5+21*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>247.52 && z_real-Z0[layer]-DZ_inWidth<=262.56)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-247.52)/pitchC4[7]))+0.5+22*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>262.56 && z_real-Z0[layer]-DZ_inWidth<=279.52)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-262.56)/pitchC4[8]))+0.5+23*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>279.52 && z_real-Z0[layer]-DZ_inWidth<=298.72)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-279.52)/pitchC4[9]))+0.5+24*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>298.72 && z_real-Z0[layer]-DZ_inWidth<=320.32)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-298.72)/pitchC4[10]))+0.5+25*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>320.32 && z_real-Z0[layer]-DZ_inWidth<=344.80)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-320.32)/pitchC4[11]))+0.5+26*32);
		 }
		 else if(z_real-Z0[layer]-DZ_inWidth>344.80 && z_real-Z0[layer]-DZ_inWidth<=372.32)
		 {
		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-344.80)/pitchC4[12]))+0.5+27*32);
		 } 
	   }
	 else if(layer==3) 
	   { // this for "C5" layers, i.e. measuring z
	     z_real = (double) (G4RandGauss::shoot(z,sigma_td));
	     ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth)/pitchC5[0]))+0.5);
	   }
  	 else if(layer==5) 
  	   { // this for "C6" layers, i.e. measuring z
  	  	 z_real = (double) (G4RandGauss::shoot(z,sigma_td));
  		 if(z_real-Z0[layer]-DZ_inWidth>0 && z_real-Z0[layer]-DZ_inWidth<=17.28)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth)/pitchC6[0]))+0.5);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>17.28 && z_real-Z0[layer]-DZ_inWidth<=32.64)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-17.28)/pitchC6[1]))+0.5+32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>32.64 && z_real-Z0[layer]-DZ_inWidth<=46.4)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-32.64)/pitchC6[2]))+0.5+2*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>46.4 && z_real-Z0[layer]-DZ_inWidth<=58.88)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-46.4)/pitchC6[3]))+0.5+3*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>58.88 && z_real-Z0[layer]-DZ_inWidth<=291.2)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-58.88)/pitchC6[4]))+0.5+4*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>291.2 && z_real-Z0[layer]-DZ_inWidth<=312.96)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-291.2)/pitchC6[5]))+0.5+26*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>312.96 && z_real-Z0[layer]-DZ_inWidth<=325.12)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-312.96)/pitchC6[6]))+0.5+28*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>325.12 && z_real-Z0[layer]-DZ_inWidth<=338.24)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-325.12)/pitchC6[7]))+0.5+29*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>338.24 && z_real-Z0[layer]-DZ_inWidth<=352.64)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-338.24)/pitchC6[8]))+0.5+30*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>352.64 && z_real-Z0[layer]-DZ_inWidth<=368.32)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-352.64)/pitchC6[9]))+0.5+31*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>368.32 && z_real-Z0[layer]-DZ_inWidth<=385.28)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-368.32)/pitchC6[10]))+0.5+32*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>385.28 && z_real-Z0[layer]-DZ_inWidth<=403.52)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-385.28)/pitchC6[11]))+0.5+33*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>403.52 && z_real-Z0[layer]-DZ_inWidth<=423.36)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-403.52)/pitchC6[12]))+0.5+34*32);
  		 }
  		 else if(z_real-Z0[layer]-DZ_inWidth>423.36 && z_real-Z0[layer]-DZ_inWidth<=444.8)
  		 {
  		 	ClosestStrip = (int) (floor(((z_real-Z0[layer]-DZ_inWidth-423.36)/pitchC6[13]))+0.5+35*32);
  		 }	 
  	   }
	   
	 else if(layer==1)
	   { // this for "Z4" layers, i.e. measuring phi
	     phi_real = phi + ((G4RandGauss::shoot(0,sigma_td))/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer]; // the sign of the 2nd term (Lorentz angle) should be a "-" as the B field is along +z and the MM are convex (and electrons are negatively charged)
	     ClosestStrip = (int) (floor(((R[layer]/pitchZ4)*(phi_real-phiij+Pi/Nsector - (Inactivtheta[layer]/2.)*Pi/180. - DZ_inLength/R[layer]))+0.5));
	   }
  	 else if(layer==2)
  	   { // this for "Z5" layers, i.e. measuring phi
  	     phi_real = phi + ((G4RandGauss::shoot(0,sigma_td))/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer]; // the sign of the 2nd term (Lorentz angle) should be a "-" as the B field is along +z and the MM are convex (and electrons are negatively charged)
  	     ClosestStrip = (int) (floor(((R[layer]/pitchZ5)*(phi_real-phiij+Pi/Nsector - (Inactivtheta[layer]/2.)*Pi/180. - DZ_inLength/R[layer]))+0.5));
  	   }
     else if(layer==4)
       { // this for "Z6" layers, i.e. measuring phi
    	phi_real = phi + ((G4RandGauss::shoot(0,sigma_td))/cos(theta_L)-(sqrt(x*x+y*y)-R[layer]+hStrip2Det)*tan(theta_L))/R[layer]; // the sign of the 2nd term (Lorentz angle) should be a "-" as the B field is along +z and the MM are convex (and electrons are negatively charged)
    	ClosestStrip = (int) (floor(((R[layer]/pitchZ6)*(phi_real-phiij+Pi/Nsector - (Inactivtheta[layer]/2.)*Pi/180. - DZ_inLength/R[layer]))+0.5));
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

