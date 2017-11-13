// gemc headers
#include "fmt_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>

vector<double> fmt_strip::FindStrip(int layer, int sector, double x, double y, double z, double Edep, fmtConstants fmtc)
{
  
	// the return vector is always in pairs.
	// The first number is the ID,
	// the second number is the sharing percentage
	vector<double> strip_id;
	// number of electrons (Nt)
	Nel = (int) (1e6*Edep/fmtc.w_i);
	sigma_td = fmtc.SigmaDrift*(z-fmtc.Z0[layer]); // expression without Lorentz angle

	//Lorentz Angle correction
	x=x+(z-fmtc.Z0[layer])*tan(fmtc.ThetaL)*cos(fmtc.Theta_Ls);
	y=y+(z-fmtc.Z0[layer])*tan(fmtc.ThetaL)*sin(fmtc.Theta_Ls);

	if (fmtc.HV_DRIFT[layer]==0||(fmtc.HV_STRIPS_IN[layer]==0&&sqrt(x*x+y*y)<fmtc.R_IR)||(fmtc.HV_STRIPS_OUT[layer]==0&&sqrt(x*x+y*y)>=fmtc.R_IR&&sqrt(x*x+y*y)<fmtc.R_max)) Nel=0;
	
	int ClosestStrip=0;
	int clust_size=fmtc.nb_sigma*sigma_td/fmtc.pitch+1; //+1 might be a little bit conservative
	if(Nel>0&&sqrt(x*x+y*y)<fmtc.R_max && sqrt(x*x+y*y)>fmtc.R_min)
	{
	
	  x_real = x*cos(fmtc.alpha[layer])+y*sin(fmtc.alpha[layer]);
	  y_real = y*cos(fmtc.alpha[layer])-x*sin(fmtc.alpha[layer]);
	  
	  if(y_real>-fmtc.y_central && y_real < fmtc.y_central){ 
	    if (x_real<=0) ClosestStrip = (int) (floor((fmtc.y_central-y_real)/fmtc.pitch)+1);
	    if (x_real>0) ClosestStrip = (int) (floor((y_real+fmtc.y_central)/fmtc.pitch)+1) + fmtc.N_halfstr+fmtc.N_sidestr;
	  }
	  else if(y_real <= -fmtc.y_central && y_real > -fmtc.R_max){ 
	    ClosestStrip = (int) (floor((fmtc.y_central-y_real)/fmtc.pitch)+1); 
	  }
	  else if(y_real >= fmtc.y_central && y_real < fmtc.R_max){ 
	    ClosestStrip = (int) (floor((y_real+fmtc.y_central)/fmtc.pitch)+1) + fmtc.N_halfstr+fmtc.N_sidestr;  
	  }
	 
	  int strip_num=ClosestStrip;// To look around closeststrip
	  double weight=Weight_td(ClosestStrip, x_real, y_real, z, fmtc);
	  if(ClosestStrip>=1 && ClosestStrip<=fmtc.N_str)
	  { // store the closest and check neighborhood of the closest (~)
	   
	    for (int clus=0; clus<clust_size+1; clus++){
	      //Look at the next strip
	      strip_num=ClosestStrip+clus;
	      if ((strip_num<fmtc.N_halfstr+fmtc.N_sidestr+1)||(ClosestStrip>=(fmtc.N_halfstr+fmtc.N_sidestr+1)&&strip_num>=(fmtc.N_halfstr+fmtc.N_sidestr+1)&&strip_num<fmtc.N_str+1)){//Check if the strip exist or make sense to look at it
		weight=Weight_td(strip_num, x_real, y_real, z, fmtc);
		strip_id.push_back(strip_num);
		strip_id.push_back(weight);
		
		if (abs(strip_y)<fmtc.y_central){//If in central part, then check the top/bottom strip systematically
		  strip_num=2*fmtc.N_halfstr+fmtc.N_sidestr+1-strip_num;
		  weight=Weight_td(strip_num, x_real, y_real, z, fmtc);
		  strip_id.push_back(strip_num);
		  strip_id.push_back(weight);
		}
	      }

	      //Look at the previous strip
	      if (clus!=0){ //Avoid double counting
		strip_num=ClosestStrip-clus;
		if (strip_num<1||(strip_num<=fmtc.N_halfstr+fmtc.N_sidestr&&ClosestStrip>fmtc.N_halfstr+fmtc.N_sidestr)) strip_num=2*fmtc.N_halfstr+fmtc.N_sidestr+1-strip_num; //Deal with the discontinuity between strip 1 and 833, or 321-513.
		weight=Weight_td(strip_num, x_real, y_real, z, fmtc);
		strip_id.push_back(strip_num);
		strip_id.push_back(weight);
		
		if (abs(strip_y)<fmtc.y_central){//If in central part, then check the top/bottom strip systematically
		  strip_num=2*fmtc.N_halfstr+fmtc.N_sidestr+1-strip_num;
		  weight=Weight_td(strip_num, x_real, y_real, z, fmtc);
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
	    for (unsigned int i=0;i<strip_id.size()/2;i++){
	      if (Nel_left==0||renorm==1){
		strip_id.at(2*i+1)=0;
	      }
	      if (renorm!=1&&Nel_left!=0){
		weight_this_strip=strip_id.at(2*i+1)/(1-renorm);
		Nel_this_strip=GetBinomial(Nel_left,weight_this_strip);
		if (Nel_this_strip==-1) cout<<Nel_left<<" -1  "<<weight_this_strip<<endl;
		renorm+=strip_id.at(2*i+1);
		strip_id.at(2*i+1)=Nel_this_strip;
		Nel_left-=Nel_this_strip;
	      }
	     
	    }
	    if (Nel_left<0) cout<<"Warning in FMT hitprocess!!!!!!!!"<<endl;
	  }
	  else
	  {// Warning - something is wrong
	    //cout<<"WARNING!!!!!! Something is wrong in FMT strip finder..... "<<ClosestStrip<<" "<<x_real<<" "<<y_real<<endl;
	    strip_id.push_back(-1);
	    strip_id.push_back(1);
	  }
	}
	if (strip_id.size()==0)
	{ // Nel=0, consider the Edep is 0
		strip_id.push_back(-1);
		strip_id.push_back(1);
	}

	return strip_id;
}

void fmt_strip::Carac_strip(int strip, fmtConstants fmtc){
  if (strip<=fmtc.N_str/2){
    strip_y=fmtc.y_central-(strip-0.5)*fmtc.pitch;
    }
    else{
      strip_y=-fmtc.y_central+(strip-fmtc.N_str/2-0.5)*fmtc.pitch;
    }
    
    //Give the strip length
    strip_length=2*fmtc.R_max*sin(acos(fabs(strip_y)/fmtc.R_max));
    if ((strip>fmtc.N_halfstr&&strip<fmtc.N_halfstr+fmtc.N_sidestr+1)||(strip<fmtc.N_str+1&&strip>2*fmtc.N_halfstr+fmtc.N_sidestr)) {
      //strip_length=2*fmtc.R_max*sin(acos(fabs(strip_y)/fmtc.R_max));
      strip_x=0;
    }	
    else{
      //if (fabs(strip_length)/fmtc.R_min<1){ 
	if (fabs(strip_y)/fmtc.R_min<1){
	strip_length=fmtc.R_max*sin(acos(fabs(strip_y)/fmtc.R_max))-fmtc.R_min*sin(acos(fabs(strip_y)/fmtc.R_min));
	if (strip<=fmtc.N_str/2) 
	  strip_x=-strip_length/2.-fmtc.R_min*sin(acos(fabs(strip_y)/fmtc.R_min));
	else
	  strip_x=strip_length/2.+fmtc.R_min*sin(acos(fabs(strip_y)/fmtc.R_min));
      }
      else{ 
	strip_length=fmtc.R_max*sin(acos(fabs(strip_y)/fmtc.R_max));
	if (strip<=fmtc.N_str/2) 
	  strip_x=-strip_length/2.;
	else
	  strip_x=strip_length/2.;
      }
    }
}


double fmt_strip::Weight_td(int strip, double x, double y, double z, fmtConstants fmtc){
  Carac_strip(strip,fmtc);
  double wght=(erf((strip_y+fmtc.pitch/2.-y)/sigma_td/sqrt(2))-erf((strip_y-fmtc.pitch/2.-y)/sigma_td/sqrt(2)))*(erf((strip_x+strip_length/2.-x)/sigma_td/sqrt(2))-erf((strip_x-strip_length/2.-x)/sigma_td/sqrt(2)))/2./2.;
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






