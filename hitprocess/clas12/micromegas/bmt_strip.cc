// gemc headers
#include "bmt_strip.h"
#include "Randomize.hh"

// c++ headers
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>

// Veronique Ziegler (Dec. 3 2015)
// Note: this method only contains the constants for the third micromegas region,
// the constants for the other 2 are missing...
// M. Ungaro (Jan 26 2016)
// Updated to read constants at the beginning of each run
// M. Defurne (May 5 2017)
// Removed the loop over the number of electrons


// the routine to find the strip
vector<double> bmt_strip::FindStrip(int layer, int sector, G4ThreeVector xyz, double Edep, bmtConstants bmtc)
{
	double x = xyz.x()/mm;
	double y = xyz.y()/mm;
	double z = xyz.z()/mm;

	int Nel = (int) (1e6*Edep/bmtc.w_i);
	// the return vector is always in pairs the first index is the strip number, the second is the Edep on the strip
	if (bmtc.HV_DRIFT[layer-1][sector-1]==0||bmtc.HV_STRIPS[layer-1][sector-1]==0) Nel=0;
	vector<double> strip_id;

	double Delta_drift = sqrt(x*x+y*y) - bmtc.RADIUS[layer-1];

	double phi = atan2(y,x) + Delta_drift*tan(bmtc.ThetaL)*cos(bmtc.Theta_Ls_Z)/bmtc.RADIUS[layer-1]; // Already apply the Lorentz Angle to find the ClosestStrip
	z=z + Delta_drift * tan(bmtc.ThetaL) * cos(bmtc.Theta_Ls_C); //Not sure useful, but take into account LorentzAngle deviation if
	int sector_bis=isInSector(layer,atan2(y,x),bmtc);
	int strip_num = getClosestStrip(layer, sector_bis, phi, z, bmtc);
	sigma = getSigma(layer, x, y, bmtc);
	sigma_phi = sigma/bmtc.RADIUS[layer-1];
	int cluster_size=0;
	if (bmtc.AXIS[layer-1]==0) cluster_size=bmtc.nb_sigma*sigma/bmtc.PITCH[layer-1][5]+1; //Compare to smallest pitch in C
	if (bmtc.AXIS[layer-1]==1) cluster_size=bmtc.nb_sigma*sigma_phi/bmtc.PITCH[layer-1][0]+1;

	double weight=0;
	if (strip_num>=1&&strip_num<=bmtc.NSTRIPS[layer-1]) weight=Weight_td(layer, sector_bis, strip_num, phi, z, bmtc);

	if(Nel>0&&weight>0) // if the track deposited energy is greater than the assumed ionization potential digitize
	  {
	    strip_id.push_back(strip_num);
	    strip_id.push_back(weight);
	   
	    // if the strip is found (i.e. the hit is within acceptance
	    for(int istrip=1;istrip< cluster_size+1;istrip++)
	      {
		//Check the strip after the closest strip
		if (strip_num+istrip<=bmtc.NSTRIPS[layer-1]) {
		  weight=Weight_td(layer, sector_bis, strip_num+istrip, phi, z, bmtc);
		  if (weight>0){
		    strip_id.push_back(strip_num+istrip);
		    strip_id.push_back(weight);
		    
		  }
		}
		//Check the strip before the closest strip
		if (strip_num-istrip>=1) {
		  weight=Weight_td(layer, sector_bis, strip_num-istrip, phi, z, bmtc);
		  if (weight>0){
		    strip_id.push_back(strip_num-istrip);
		    strip_id.push_back(weight);
		   }
		}
		
	      }

	    //We have computed the weight with a gaussian distribution. But considering the few electrons, it makes no sense.
	    int Nel_left=Nel;
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
		//cout<<layer<<" "<<sector_bis<<" "<<strip_id.at(2*i)<<"  "<<strip_id.at(2*i+1)<<endl;
	      }
	    }

	    //Final check... Energy must be conserved unfortunately
	    if (Nel_left<0) cout<<"Warning in BMT hitprocess!!!!!!!!"<<endl;
	  }
	if (strip_id.size()==0)
	  { // Nel=0, consider the Edep is 0
	    strip_id.push_back(-1);
	    strip_id.push_back(-1);
	  }
       
	return strip_id;
}


//
// param layer
// param x x-coordinate of the hit in the lab frame
// param y y-coordinate of the hit in the lab frame
// return the sigma in the azimuth direction taking the Lorentz angle into account
//
double bmt_strip::getSigma(int layer, double x, double y,  bmtConstants bmtc)
{ // sigma for Z-detectors

	double sigma = bmtc.SigmaDrift*(sqrt(x*x+y*y) - bmtc.RADIUS[layer-1])/cos(bmtc.ThetaL);

	return sigma;

}

int bmt_strip::getClosestStrip(int layer, int sector, double angle, double z, bmtConstants bmtc){
  double var=0;
  double var_min=0;
  double var_max=0; //var=z if it is a C detector, var=angle if Z detector
  int group=0;
  double strip_offset=0; //For C because of the variable pitch
  int ClosestStrip=-1;
 
  if(angle<0) angle+=2*pi; // from 0 to 2Pi
 
  //To deal with the sector covering 0 degree
  double angle_i = bmtc.EDGE1[layer-1][sector]; // first angular boundary of sector
  double angle_f = bmtc.EDGE2[layer-1][sector]; // second angular boundary of sector
 
  if (angle_f<angle_i) { //We are on sector covering 0 degree
    angle_f+= 2*pi;
    if (angle<angle_i) angle+=2*pi; //If angle smaller than angle_i, then we have passed 0 degree
  }
 
  if (angle>angle_i&&angle<angle_f&&z<bmtc.ZMAX[layer-1]&&z>bmtc.ZMIN[layer-1]){
    if (bmtc.AXIS[layer-1]==0){//Then it is a C detector
      var=z;
      var_min=bmtc.ZMIN[layer-1];
      var_max=var_min+bmtc.GROUP[layer-1][0]*bmtc.PITCH[layer-1][0];
    }
    if (bmtc.AXIS[layer-1]==1){//Then it is a Z detector
      var=angle;
      var_min=angle_i;
      var_max=angle_f;
    }
    while (var>var_max){
      strip_offset+=bmtc.GROUP[layer-1][group];
      var_min+=bmtc.GROUP[layer-1][group]*bmtc.PITCH[layer-1][group];
      var_max+=bmtc.GROUP[layer-1][group+1]*bmtc.PITCH[layer-1][group+1];
      group++;
    }
    ClosestStrip=strip_offset+floor((var-var_min)/bmtc.PITCH[layer-1][group])+1;
  }
  return ClosestStrip;
}

// Return the group of equally separated strips in which the strip is
int bmt_strip::getStripGroup(int layer, int strip, bmtConstants bmtc){
  int group=0; //Z is always one group
  int total_strip=bmtc.GROUP[layer-1][group];
  while (strip>total_strip){
    group++;
    total_strip+=bmtc.GROUP[layer-1][group];
  }
  return group;
}

// param layer the hit layer
// param strip the hit strip
// return the z position in mm for the C-detectors
//
// WARNING: This routine is not used?
//
double bmt_strip::GetStripInfo(int layer, int sector, int strip, bmtConstants bmtc)
{
	int num_strip = strip - 1;     			// index of the strip (starts at 0)
	double var=0.;
	if (bmtc.AXIS[layer-1]==0) var=bmtc.ZMIN[layer-1]; //C detector so we look at Z
	if (bmtc.AXIS[layer-1]==1) var=bmtc.EDGE1[layer-1][sector]; //Z detector so we look at phi

	int group=0;
	int limit = bmtc.GROUP[layer-1][group];

	if (num_strip>0){
	  for (int j=1;j<num_strip;j++)
	  {
	    var += bmtc.PITCH[layer-1][group];
		if (j>=limit)
		{ //test if we change the width
			group++;
			limit += bmtc.GROUP[layer-1][group];
		}
		var += bmtc.PITCH[layer-1][group];
	  }
	}
	var+=0.5*bmtc.PITCH[layer-1][group];
	return var; 
}


// not used yet (implemented for alternate algorithm not yet fully developed)
int bmt_strip::isInSector(int layer, double angle, bmtConstants bmtc)
{
	if(angle<0)
		angle+=2*pi; // from 0 to 2Pi
	double angle_pr = angle + 2*pi;

	double angle_i = 0; // first angular boundary init
	double angle_f = 0; // second angular boundary for detector A, B, or C init
	int num_detector = -1;

	for(int i = 0; i<3; i++) {

		angle_i = bmtc.EDGE1[layer-1][i] ;
		angle_f = bmtc.EDGE2[layer-1][i] ;
		if (angle_f<angle_i) angle_f+= 2*pi;

		if((angle>=angle_i && angle<=angle_f)||(angle_pr>=angle_i && angle_pr<=angle_f)) num_detector=i;
	}
	return num_detector;
}

double bmt_strip::Weight_td(int layer, int sector, int strip, double angle, double z, bmtConstants bmtc){
  double wght=0;
  int group=getStripGroup(layer, strip, bmtc);
  if(angle<0) angle+=2*pi; // from 0 to 2Pi
  
  //To deal with the sector covering 0 degree
  double angle_i = bmtc.EDGE1[layer-1][sector]; // first angular boundary of sector
  double angle_f = bmtc.EDGE2[layer-1][sector]; // second angular boundary of sector
  if (angle_f<angle_i) { //We are on sector covering 0 degree
    angle_f+= 2*pi;
    if (angle<angle_i) angle+=2*pi; //If angle smaller than angle_i, then we have passed 0 degree
  }
  
  if(bmtc.AXIS[layer-1]==1){ // if it is a Z-detector
    double strip_phi=angle_i+(strip-0.5)*bmtc.PITCH[layer-1][0];
    double strip_z=(bmtc.ZMIN[layer-1]+bmtc.ZMAX[layer-1])/2.;
    double strip_length=bmtc.ZMAX[layer-1]-bmtc.ZMIN[layer-1];
    wght=(erf((strip_phi+bmtc.PITCH[layer-1][0]/2.-angle)/sigma_phi/sqrt(2))-erf((strip_phi-bmtc.PITCH[layer-1][0]/2.-angle)/sigma_phi/sqrt(2)))*(erf((strip_z+strip_length/2.-z)/sigma/sqrt(2))-erf((strip_z-strip_length/2.-z)/sigma/sqrt(2)))/2./2.;
  }
  if(bmtc.AXIS[layer-1]==0){ // if it is a C-detector
    double strip_phi=(angle_i+angle_f)/2.;
    double strip_z=bmtc.ZMIN[layer-1];
    int strip_offset=0;
    for (int i=0;i<group;i++){
      strip_z+=bmtc.GROUP[layer-1][i]*bmtc.PITCH[layer-1][i];
      strip_offset+=bmtc.GROUP[layer-1][i];
    }
    strip_z+=(strip-strip_offset-0.5)*bmtc.PITCH[layer-1][group];
    double strip_length=bmtc.EDGE2[layer-1]-bmtc.EDGE1[layer-1];
    wght=(erf((strip_z+bmtc.PITCH[layer-1][group]/2.-z)/sigma/sqrt(2))-erf((strip_z-bmtc.PITCH[layer-1][group]/2.-z)/sigma/sqrt(2)))*(erf((strip_phi+strip_length/2.-angle)/sigma_phi/sqrt(2))-erf((strip_phi-strip_length/2.-angle)/sigma_phi/sqrt(2)))/2./2.;
  }
  if (wght<0) wght=-wght;
  return wght;
}

double bmt_strip::GetBinomial(double n, double p){
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





