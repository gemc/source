// gemc headers
#include "recoil_strip.h"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

vector<recoil_strip_found> recoil_strip::FindStrip(G4ThreeVector xyz , double Edep, recoilConstants recoilc, double time)
{	
	vector<recoil_strip_found> strip_found;
	vector<recoil_strip_found> strip_found_temp;
	recoil_strip_found ClosestStrip;
	double time_strip =0; // gauss(time_gemc + time_dz + time_redout, sigma_dt);
		
	int N_el = 1e6*Edep/recoilc.w_i;

	if (N_el ==0){
		ClosestStrip.numberID = -15000;
		ClosestStrip.weight = 1;
		ClosestStrip.time = -1;
		strip_found.push_back(ClosestStrip);
		return strip_found;
	}
	
	N_el = G4Poisson(N_el*recoilc.gain);
	
	// strip reference frame
		
       	double x_real = xyz.x()*cos(M_PI*recoilc.get_stereo_angle()/180) + xyz.y()*sin(M_PI*recoilc.get_stereo_angle()/180);
	double y_real = xyz.y()*cos(M_PI*recoilc.get_stereo_angle()/180) - xyz.x()*sin(M_PI*recoilc.get_stereo_angle()/180); 
	double z_real = xyz.z(); 
	
	double time_dz = fabs(-recoilc.Zhalf/cm + z_real/cm)/recoilc.v_drift ;
	int ClosestStrip_ID;

	ClosestStrip_ID=round((y_real)/recoilc.get_strip_pitch());

	double weight=Weight_td(ClosestStrip_ID, x_real, y_real, z_real, recoilc);

	double strip_length_toReadout =cal_length(strip_endpoint1, xyz);
	double time_toReadout = strip_length_toReadout/recoilc.v_eff_readout;
	if(recoilc.v_eff_readout ==0) time_toReadout=0;
	time_strip = time + time_dz + time_toReadout + G4RandGauss::shoot(0., recoilc.sigma_time);
	
	ClosestStrip.numberID = ClosestStrip_ID;
	ClosestStrip.weight = weight;
	ClosestStrip.time = time_strip;
	strip_found_temp.push_back(ClosestStrip);
	
	//To look around closest strip
	recoil_strip_found NextStrip;
	
	int strip_num=0;
	double weight_next=1;
	double weight_previous=1;
	int clus =1;

	if(recoilc.get_strip_kind()=="strip_u") {number_of_strip=round(2*recoilc.Yhalf/recoilc.get_strip_pitch());}
	if(recoilc.get_strip_kind()=="strip_v") {number_of_strip=round(2*recoilc.Xhalf/recoilc.get_strip_pitch());}
	
	while((weight_next>=0. || weight_previous>=0.)){
	//Look at the next strip
		strip_num = ClosestStrip_ID + clus;
		weight_next = Weight_td(strip_num, x_real, y_real, z_real, recoilc);
		if(weight_next!=-1){
		  NextStrip.numberID = strip_num;
		  NextStrip.weight = weight_next;
		  strip_length_toReadout =cal_length(strip_endpoint1, xyz);
		  time_toReadout = strip_length_toReadout/recoilc.v_eff_readout;
		  if(recoilc.v_eff_readout ==0) time_toReadout=0;
		  NextStrip.time = time + time_dz + time_toReadout + G4RandGauss::shoot(0., recoilc.sigma_time);
		  strip_found_temp.push_back(NextStrip);
		}		
		//Look at the previous strip
		strip_num = ClosestStrip_ID - clus;
		weight_previous = Weight_td(strip_num, x_real, y_real, z_real, recoilc);
		if(weight_previous!=-1){
		  NextStrip.numberID = strip_num;
		  NextStrip.weight = weight_previous;
		  strip_length_toReadout =cal_length(strip_endpoint1, xyz);
		  time_toReadout = strip_length_toReadout/recoilc.v_eff_readout;
		  if(recoilc.v_eff_readout ==0) time_toReadout=0;
		  NextStrip.time = time + time_dz + time_toReadout + G4RandGauss::shoot(0., recoilc.sigma_time);
		  strip_found_temp.push_back(NextStrip);
		}
		clus++;
	}
	
	/* New strip ID numeration: 1.... Number_of_strips. Number of involved strips is the size of the vector strip_found_temp  */
	
	auto max = std::max_element( strip_found_temp.begin(), strip_found_temp.end(),
				     []( const recoil_strip_found &a, const recoil_strip_found &b )
				     {
				       return a.numberID < b.numberID;
				     } );
	
	auto min = std::min_element( strip_found_temp.begin(), strip_found_temp.end(),
				     []( const recoil_strip_found &a, const recoil_strip_found &b )
				     {
				       return a.numberID < b.numberID;
				     } );
	
       	auto avg = round((max->numberID + min->numberID)/2);

       	for (unsigned int i=0; i<strip_found_temp.size();i++){
	  strip_found_temp.at(i).numberID = strip_found_temp.at(i).numberID - avg + strip_found_temp.size()/2 ;
       	}
	
	double Nel_left=N_el;
	double renorm=0;
	double weight_this_strip;
	int Nel_this_strip=0;
	
	for (unsigned int i=0;i<strip_found_temp.size();i++){
		if (Nel_left==0||renorm==1){
			strip_found_temp.at(i).weight=0;
		}
		if (renorm!=1&&Nel_left!=0){
			weight_this_strip=strip_found_temp.at(i).weight/(1-renorm);
			Nel_this_strip=GetBinomial(Nel_left,weight_this_strip);
			renorm+=strip_found_temp.at(i).weight;
			strip_found_temp.at(i).weight=Nel_this_strip;
			Nel_left-=Nel_this_strip;
		}
	}
	
	for (unsigned int i=0;i<strip_found_temp.size();i++){
		if(strip_found_temp.at(i).weight>0) strip_found.push_back(strip_found_temp[i]);
	}
	
	if (strip_found.size() ==0){
		ClosestStrip.numberID = -15000;
		ClosestStrip.weight = 1;
		ClosestStrip.time = -1;
		strip_found.push_back(ClosestStrip);		
	}
	/*for (unsigned int i=0;i<strip_found.size();i++){
	  	  if(strip_found.at(i).numberID>number_of_strip||(strip_found.at(i).numberID<0 && strip_found.at(i).numberID!=-15000)) cout<<"ATTENTION! "<<"x "<<x_real<<" y "<<y_real<<" "<<recoilc.get_strip_kind()<<" numberID "<<strip_found.at(i).numberID<<" c "<<strip_found.at(i).numberID*recoilc.get_strip_pitch()<<" 2*base_x "<<2*recoilc.Xhalf<<" weight "<<strip_found.at(i).weight<<" max "<<number_of_strip<<endl;
		  }*/
	return strip_found;
	
}

double recoil_strip::Weight_td(int strip, double x, double y, double z, recoilConstants recoilc){
	double wght;
	if(Build_strip(strip, recoilc)){
	 wght=(erf((strip_y+recoilc.get_strip_width(strip)/2.-y)/recoilc.sigma_td/sqrt(2))-erf((strip_y-recoilc.get_strip_width(strip)/2.-y)/recoilc.sigma_td/sqrt(2)))*(erf((strip_x+strip_length/2.-x)/recoilc.sigma_td/sqrt(2))-erf((strip_x-strip_length/2.-x)/recoilc.sigma_td/sqrt(2)))/2./2.;
	 if (wght<0) wght=-wght;
	}else{
		wght =-1;
	}
	return wght;
}

bool recoil_strip::Build_strip(int strip, recoilConstants recoilc){

        double c = strip*recoilc.get_strip_pitch();

   // Trapezoid coordinates
	G4ThreeVector A = {-recoilc.Xhalf, -recoilc.Yhalf, recoilc.Zhalf};
	G4ThreeVector B =  {recoilc.Xhalf, -recoilc.Yhalf, recoilc.Zhalf};
	G4ThreeVector C = {-recoilc.Xhalf, recoilc.Yhalf, recoilc.Zhalf};
	G4ThreeVector D =  {recoilc.Xhalf, recoilc.Yhalf, recoilc.Zhalf};
	
	// C-------------D //
	// --------------- //
	// --------------- //
	// A-------------B //
	// Intersection points between strip straight line and Trapezoid straight lines

	G4ThreeVector AB_strip={0.,0.,0.};
	G4ThreeVector BD_strip={0.,0.,0.};
	G4ThreeVector CD_strip={0.,0.,0.};
	G4ThreeVector AC_strip={0.,0.,0.};
	
	if(recoilc.get_strip_kind()=="strip_v" && fabs(c)<=recoilc.Xhalf) {
	  AB_strip = {c,-recoilc.Yhalf,recoilc.Zhalf};
	  CD_strip = {c,recoilc.Yhalf,recoilc.Zhalf};
	}
	else if(recoilc.get_strip_kind()=="strip_u" && fabs(c)<=recoilc.Yhalf) {
	  BD_strip = {recoilc.Xhalf,c,recoilc.Zhalf};
	  AC_strip = {-recoilc.Xhalf,c,recoilc.Zhalf};
	}
	else return false;	
	
	// geometrical characteristic
	double length_strip=0;
	G4ThreeVector first_point;
	G4ThreeVector second_point;
	
	// check if the intersection point is on the segment defined by two points (i.e A and B)
	
	 if(recoilc.get_strip_kind()=="strip_v"){
	   first_point=AB_strip;
	   second_point=CD_strip;
	}	
     	
	 if(recoilc.get_strip_kind()=="strip_u"){
           first_point=AC_strip;
	   second_point=BD_strip;		
	}
	
	length_strip = cal_length(first_point, second_point);
	strip_length = length_strip;
	if(recoilc.get_strip_kind()=="strip_u") strip_length= 2*recoilc.Xhalf;
	if(recoilc.get_strip_kind()=="strip_v") strip_length= 2*recoilc.Yhalf;

	strip_endpoint1 = first_point;
	strip_endpoint2 = second_point;

	G4ThreeVector strip_endpoint1_stripFrame = change_of_coordinates(strip_endpoint1, recoilc);
	G4ThreeVector strip_endpoint2_stripFrame = change_of_coordinates(strip_endpoint2, recoilc);
	
	strip_y = strip_endpoint1_stripFrame.y();
	strip_x = (strip_endpoint1_stripFrame.x() + strip_endpoint2_stripFrame.x())/2;

	return true;

}

G4ThreeVector recoil_strip::change_of_coordinates( G4ThreeVector A, recoilConstants recoilc){
	
	G4ThreeVector XYZ;
	XYZ.setX(A.x()*cos(M_PI*recoilc.get_stereo_angle()/180) + A.y()*sin(M_PI*recoilc.get_stereo_angle()/180));
        XYZ.setY(A.y()*cos(M_PI*recoilc.get_stereo_angle()/180) - A.x()*sin(M_PI*recoilc.get_stereo_angle()/180));
	XYZ.setZ(A.z());
	
	return XYZ;
}

double recoil_strip::cal_length(G4ThreeVector A, G4ThreeVector B){
	double length=0;
	length = sqrt(pow((A.x()-B.x()),2) + pow((A.y()-B.y()),2)) ;
	//cout<<"length in cal_length: "<<length<<endl;
	return length;
}

double recoil_strip::GetBinomial(double n, double p){
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
