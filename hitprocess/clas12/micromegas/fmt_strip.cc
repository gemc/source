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
	intersuperlayer = 11.9;      // distance between 2 superlayers // modified on 7/27/2015 to match new geometry (Frederic Georges)
	//  pitch           = 0.500;     // pitch of the strips
	pitch           = 0.525;     // pitch of the strips // modified on 7/27/2015 to match new geometry (Frederic Georges)
	hDrift          = 5.0;
	//hStrip2Det      = hDrift/2.;
	sigma_td_max    = 0.01; // very small transverse diffusion because of B field
	w_i             = 25.0;

	R_min = 42.575;   // modified on 7/27/2015 to match new geometry (Frederic Georges)
  // outer radius of strip part // modified on 7/27/2015 to match new geometry (Frederic Georges)
	//  Z_1stlayer = 295.-0.5-0.1-0.015-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap)
	//Z_1stlayer = 305.25-0.5-0.1-0.015-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap) // modified on 5/20/2015 to match new geometry (R.De Vita)
	//Z_1stlayer = 305.250-0.005-0.200-2.000-0.200-0.005-0.050-0.020-0.128-0.030-2.500; // z position of the 1st layer (middle of the drift gap) // modified on 7/27/2015 to match new geometry (Frederic Georges)
	Z_1stlayer = 300.3;
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
	Z0.push_back(Z0[2]+intersuperlayer+2.0); //modified on 2/22/2017 to match new geometry (Maxime Defurne)
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
	N_halfstr = 320; //number of bottom strips in the central part
	N_sidestr = (N_str-2*N_halfstr)/2; //number of strips one side
	y_central = N_halfstr*pitch/2.; // Y-limit of the central part
	nb_sigma=4;

	R_max = pitch*(N_halfstr+2*N_sidestr)/2.;
	
}

vector<double> fmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep)
{
	// the return vector is always in pairs.
	// The first number is the ID,
	// the second number is the sharing percentage
	vector<double> strip_id;
	// number of electrons (Nt)
	Nel = (int) (1e6*Edep/w_i);
	sigma_td = sigma_td_max* sqrt((z-Z0[layer])/hDrift); // expression without Lorentz angle
	//if(z>Z0[layer]+7.697 || z<Z0[layer]+2.697) cout << "Warning! z position of the FMT hit is not in the sensitive volume: " << z <<" for bottom of gas at "<<Z0[layer]<< endl;

	
	int ClosestStrip=0;
	int clust_size=nb_sigma*sigma_td/pitch+1; //+1 might be a little bit conservative
	if(Nel>0&&sqrt(x*x+y*y)<R_max && sqrt(x*x+y*y)>R_min)
	{
	
	  x_real = x*cos(alpha[layer])+y*sin(alpha[layer]);
	  y_real = y*cos(alpha[layer])-x*sin(alpha[layer]);
	  
	  if(y_real>-y_central && y_real < y_central){ 
	    if (x_real<=0) ClosestStrip = (int) (floor((y_central-y_real)/pitch)+1);
	    if (x_real>0) ClosestStrip = (int) (floor((y_real+y_central)/pitch)+1) + N_halfstr+N_sidestr;
	  }
	  else if(y_real <= -y_central && y_real > -R_max){ 
	    ClosestStrip = (int) (floor((y_central-y_real)/pitch)+1); 
	  }
	  else if(y_real >= y_central && y_real < R_max){ 
	    ClosestStrip = (int) (floor((y_real+y_central)/pitch)+1) + N_halfstr+N_sidestr;  
	  }
	 
	  int strip_num=ClosestStrip;// To look around closeststrip
	  double weight=Weight_td(ClosestStrip, x_real, y_real, z);
	  if(ClosestStrip>=1 && ClosestStrip<=N_str)
	  { // store the closest and check neighborhood of the closest (~)
	   
	    for (int clus=0; clus<clust_size+1; clus++){
	      //Look at the next strip
	      strip_num=ClosestStrip+clus;
	      if ((strip_num<N_halfstr+N_sidestr+1)||(ClosestStrip>=(N_halfstr+N_sidestr+1)&&strip_num>=(N_halfstr+N_sidestr+1)&&strip_num<N_str+1)){//Check if the strip exist or make sense to look at it
		weight=Weight_td(strip_num, x_real, y_real, z);
		strip_id.push_back(strip_num);
		strip_id.push_back(weight);
		
		if (abs(strip_y)<y_central){//If in central part, then check the top/bottom strip systematically
		  strip_num=2*N_halfstr+N_sidestr+1-strip_num;
		  weight=Weight_td(strip_num, x_real, y_real, z);
		  strip_id.push_back(strip_num);
		  strip_id.push_back(weight);
		}
	      }

	      //Look at the previous strip
	      if (clus!=0){ //Avoid double counting
		strip_num=ClosestStrip-clus;
		if (strip_num<1||(strip_num<=N_halfstr+N_sidestr&&ClosestStrip>N_halfstr+N_sidestr)) strip_num=2*N_halfstr+N_sidestr+1-strip_num; //Deal with the discontinuity between strip 1 and 833, or 321-513.
		weight=Weight_td(strip_num, x_real, y_real, z);
		strip_id.push_back(strip_num);
		strip_id.push_back(weight);
		
		if (abs(strip_y)<y_central){//If in central part, then check the top/bottom strip systematically
		  strip_num=2*N_halfstr+N_sidestr+1-strip_num;
		  weight=Weight_td(strip_num, x_real, y_real, z);
		  strip_id.push_back(strip_num);
		  strip_id.push_back(weight);
		} 
	      }
	    }
	     //We have computed the weight with a gaussian distribution. But considering the few electrons, it makes no sense.
	    double Nel_left=Nel;
	    double renorm=0;
	    double weight_this_strip;
	    int Nel_this_strip=0;
	    for (int i=0;i<strip_id.size()/2;i++){
	      if (renorm!=1&&Nel_left!=0){
		weight_this_strip=strip_id.at(2*i+1)/(1-renorm);
		Nel_this_strip=GetBinomial(Nel_left,weight_this_strip);
		if (Nel_this_strip==-1) cout<<Nel_left<<" -1  "<<weight_this_strip<<endl;
		renorm+=strip_id.at(2*i+1);
		strip_id.at(2*i+1)=Nel_this_strip;
		Nel_left-=Nel_this_strip;
	      }
	      //	cout<<"Strip number "<<strip_id.at(2*i)<<" "<<strip_id.at(2*i+1)<<endl;
	    }
	    if (Nel_left<0) cout<<"Warning in FMT hitprocess!!!!!!!!"<<endl;
	  }
	  else
	  {// Warning - something is wrong
	    cout<<"WARNING!!!!!! Something is wrong in FMT strip finder..... "<<ClosestStrip<<" "<<x_real<<" "<<y_real<<endl;
	    strip_id.push_back(-1);
	    strip_id.push_back(1);
	  }
	}
	else
	{ // Nel=0, consider the Edep is 0
		strip_id.push_back(-1);
		strip_id.push_back(1);
	}

	return strip_id;
}

void fmt_strip::Carac_strip(int strip){
  if (strip<=N_str/2){
    strip_y=y_central-(strip-0.5)*pitch;
    }
    else{
      strip_y=-y_central+(strip-N_str/2-0.5)*pitch;
    }
    
    //Give the strip length
    strip_length=2*R_max*sin(acos(fabs(strip_y)/R_max));
    if ((strip>N_halfstr&&strip<N_halfstr+N_sidestr+1)||(strip<N_str+1&&strip>2*N_halfstr+N_sidestr)) {
      //strip_length=2*R_max*sin(acos(fabs(strip_y)/R_max));
      strip_x=0;
    }	
    else{
      //if (fabs(strip_length)/R_min<1){ 
	if (fabs(strip_y)/R_min<1){
	strip_length=R_max*sin(acos(fabs(strip_y)/R_max))-R_min*sin(acos(fabs(strip_y)/R_min));
	if (strip<=N_str/2) 
	  strip_x=-strip_length/2.-R_min*sin(acos(fabs(strip_y)/R_min));
	else
	  strip_x=strip_length/2.+R_min*sin(acos(fabs(strip_y)/R_min));
      }
      else{ 
	strip_length=R_max*sin(acos(fabs(strip_y)/R_max));
	if (strip<=N_str/2) 
	  strip_x=-strip_length/2.;
	else
	  strip_x=strip_length/2.;
      }
    }
}


double fmt_strip::Weight_td(int strip, double x, double y, double z){
  Carac_strip(strip);
  double wght=(erf((strip_y+pitch/2.-y)/sigma_td/sqrt(2))-erf((strip_y-pitch/2.-y)/sigma_td/sqrt(2)))*(erf((strip_x+strip_length/2.-x)/sigma_td/sqrt(2))-erf((strip_x-strip_length/2.-x)/sigma_td/sqrt(2)))/2./2.;
  if (wght<0) wght=-wght;
  return wght;
}

double fmt_strip::GetBinomial(double n, double p){
  double answer;
  answer=CLHEP::RandBinomial::shoot(n,p);
  //Very bad method when n=0 or p close to 0 or 1... return easily -1 in these case.
  //So need to help in the limit condition
  if (answer==-1){
    answer=n;
    if (p==0) answer=0;
  }
  return answer;
}






