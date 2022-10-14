// gemc headers
#include "uRwell_strip.h"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

vector<uRwell_strip_found> uRwell_strip::FindStrip(G4ThreeVector xyz , double Edep, uRwellConstants uRwellc, double time)
{
	
	vector<uRwell_strip_found> strip_found;
	vector<uRwell_strip_found> strip_found_temp;
	uRwell_strip_found ClosestStrip;
	double time_strip =0; // gauss(time_gemc + time_dz + time_redout, sigma_dt);
	
	// int N_strip = Number_of_strip(uRwellc);
	
	int N_el = 1e6*Edep/uRwellc.w_i;
	
	if (N_el ==0){
		ClosestStrip.numberID = -15000;
		ClosestStrip.weight = 1;
		ClosestStrip.time = -1;
		strip_found.push_back(ClosestStrip);
		return strip_found;
	}
	
	N_el = G4Poisson(N_el*uRwellc.gain);
	
	// strip reference frame
	
	
	double x_real = xyz.x()*cos(M_PI*uRwellc.find_stereo_angle()/180) + xyz.y()*sin(M_PI*uRwellc.find_stereo_angle()/180);
	double y_real = xyz.y()*cos(M_PI*uRwellc.find_stereo_angle()/180) - xyz.x()*sin(M_PI*uRwellc.find_stereo_angle()/180);
	double z_real = xyz.z();
	
	
	double time_dz = fabs(-uRwellc.Zhalf/cm + z_real/cm)/uRwellc.v_drift ;
	
	int ClosestStrip_ID = round((y_real)/uRwellc.find_strip_pitch());
	
	double weight=Weight_td(ClosestStrip_ID, x_real, y_real, z_real, uRwellc);
	double strip_length_toReadout =cal_length(strip_endpoint1, xyz);
	double time_toReadout = strip_length_toReadout/uRwellc.v_eff_readout;
	if(uRwellc.v_eff_readout ==0) time_toReadout=0;
	time_strip = time + time_dz + time_toReadout + G4RandGauss::shoot(0., uRwellc.sigma_time);
	
	ClosestStrip.numberID = ClosestStrip_ID;
	ClosestStrip.weight = weight;
	ClosestStrip.time = time_strip;
	strip_found_temp.push_back(ClosestStrip);
	
	//To look around closest strip
	uRwell_strip_found NextStrip;
	
	int strip_num=0;
	double weight_next=1;
	double weight_previous=1;
	int clus =1;
	// int cont =0;
	
	while(weight_next>=0. || weight_previous>=0.){
		
		//Look at the next strip
		strip_num = ClosestStrip_ID + clus;
		weight_next = Weight_td(strip_num, x_real, y_real, z_real, uRwellc);
		if(weight_next!=-1){
			NextStrip.numberID = strip_num;
			NextStrip.weight = weight_next;
			strip_length_toReadout =cal_length(strip_endpoint1, xyz);
			time_toReadout = strip_length_toReadout/uRwellc.v_eff_readout;
			if(uRwellc.v_eff_readout ==0) time_toReadout=0;
			NextStrip.time = time + time_dz + time_toReadout + G4RandGauss::shoot(0., uRwellc.sigma_time);
			strip_found_temp.push_back(NextStrip);
		}
		
		//Look at the previous strip
		strip_num = ClosestStrip_ID - clus;
		weight_previous = Weight_td(strip_num, x_real, y_real, z_real, uRwellc);
		if(weight_previous!=-1){
			NextStrip.numberID = strip_num;
			NextStrip.weight = weight_previous;
			strip_length_toReadout =cal_length(strip_endpoint1, xyz);
			time_toReadout = strip_length_toReadout/uRwellc.v_eff_readout;
			if(uRwellc.v_eff_readout ==0) time_toReadout=0;
			NextStrip.time = time + time_dz + time_toReadout + G4RandGauss::shoot(0., uRwellc.sigma_time);
			strip_found_temp.push_back(NextStrip);
		}
		clus++;
		
	}
	
	
	/* New strip ID numeration: 1.... Number_of_strips. Number of strips is the size of the vector strip_found_temp  */
	
	auto max = std::max_element( strip_found_temp.begin(), strip_found_temp.end(),
										 []( const uRwell_strip_found &a, const uRwell_strip_found &b )
										 {
		return a.numberID < b.numberID;
	} );
	
	auto min = std::min_element( strip_found_temp.begin(), strip_found_temp.end(),
										 []( const uRwell_strip_found &a, const uRwell_strip_found &b )
										 {
		return a.numberID < b.numberID;
	} );
	
	int avg = (max->numberID + min->numberID)/2;
	//	auto number_of_strip = Number_of_strip(uRwellc);
	
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
			if (Nel_this_strip==-1) cout<<Nel_left<<" -1  "<<weight_this_strip<<endl;
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
	
	return strip_found;
	
}


double uRwell_strip::Weight_td(int strip, double x, double y, double z, uRwellConstants uRwellc){
	double wght;
	if(Build_strip(strip, uRwellc)){
		wght=( erf((strip_y+uRwellc.find_strip_width()/2.-y)/uRwellc.sigma_td/sqrt(2))-erf((strip_y-uRwellc.find_strip_width()/2.-y)/uRwellc.sigma_td/sqrt(2)))*
		(erf((strip_x+strip_length/2.-x)/uRwellc.sigma_td/sqrt(2))-erf((strip_x-strip_length/2.-x)/uRwellc.sigma_td/sqrt(2)))/2./2.;
		if (wght<0) wght=-wght;
	}else{
		wght =-1;
	}
	return wght;
}

bool uRwell_strip::Build_strip(int strip, uRwellConstants uRwellc ){
	
	//strip straight line -> y = mx +c;
	double m = tan(M_PI*uRwellc.find_stereo_angle()/180);
	double c = strip*uRwellc.find_strip_pitch()/cos(M_PI*uRwellc.find_stereo_angle()/180);
	
	// Trapezoid coordinates
	G4ThreeVector A = {uRwellc.Xhalf_base, -uRwellc.Yhalf, uRwellc.Zhalf};
	G4ThreeVector B =  {-uRwellc.Xhalf_base, -uRwellc.Yhalf, uRwellc.Zhalf};
	G4ThreeVector C = {uRwellc.Xhalf_Largebase, uRwellc.Yhalf, uRwellc.Zhalf};
	G4ThreeVector D =  {-uRwellc.Xhalf_Largebase, uRwellc.Yhalf, uRwellc.Zhalf};
	
	
	// C-------------D //
	//  -------------  //
	//   -----------   //
	//    A-------B   //
	// Intersection points between strip straight line and Trapezoid straight lines
	
	G4ThreeVector AB_strip = intersectionPoint(m,c,A,B);
	G4ThreeVector BD_strip = intersectionPoint(m,c,B,D);
	G4ThreeVector CD_strip = intersectionPoint(m,c,C,D);
	G4ThreeVector AC_strip = intersectionPoint(m,c,A,C);
	
	vector< G4ThreeVector> strip_points ; // intersection point between strip and the trapezoid sides;
	
	// geometrical characteristic
	double length_strip=0;
	// double lenght_strip_temp=0;
	G4ThreeVector first_point;
	G4ThreeVector second_point;
	
	// check if the intersection point is on the segment defined by two points (i.e A and B)
	
	if(uRwellc.find_strip_kind()=="strip_u"){
		if(pointOnsegment(AC_strip, A, C)) {
			first_point=AC_strip;
			
			if(pointOnsegment(BD_strip, B, D)) second_point = BD_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
			
		}else if(pointOnsegment(AB_strip, A, B)){
			first_point=AB_strip;
			
			if(pointOnsegment(BD_strip, B, D)) second_point = BD_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
			
		}else{
			return false;
		}
	}
	
	
	
	if(uRwellc.find_strip_kind()=="strip_v"){
		
		if(pointOnsegment(BD_strip, B, D)){
			first_point=BD_strip;
			
			if(pointOnsegment(AC_strip, A, C)) second_point = AC_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
		}else if (pointOnsegment(AB_strip, A, B)){
			first_point=AB_strip;
			if(pointOnsegment(AC_strip, A, C)) second_point = AC_strip;
			if(pointOnsegment(CD_strip, C, D)) second_point = CD_strip;
		}else{
			return false;
		}
		
	}
	
	length_strip = cal_length(first_point, second_point);
	
	strip_length = length_strip;
	strip_endpoint1 = first_point;
	strip_endpoint2 = second_point;
	
	G4ThreeVector strip_endpoint1_stripFrame = change_of_coordinates(strip_endpoint1, uRwellc );
	G4ThreeVector strip_endpoint2_stripFrame = change_of_coordinates(strip_endpoint2, uRwellc );
	
	strip_y = strip_endpoint1_stripFrame.y();
	strip_x = (strip_endpoint1_stripFrame.x() + strip_endpoint2_stripFrame.x())/2;
	
	return true;
}

G4ThreeVector uRwell_strip::change_of_coordinates( G4ThreeVector A, uRwellConstants uRwellc){
	
	G4ThreeVector XYZ;
	XYZ.setX(A.x()*cos(M_PI*uRwellc.find_stereo_angle()/180) + A.y()*sin(M_PI*uRwellc.find_stereo_angle()/180));
	XYZ.setY(A.y()*cos(M_PI*uRwellc.find_stereo_angle()/180) - A.x()*sin(M_PI*uRwellc.find_stereo_angle()/180));
	XYZ.setZ(A.z());
	
	return XYZ;
}

G4ThreeVector uRwell_strip::intersectionPoint(double m, double c, G4ThreeVector A, G4ThreeVector B){
	
	G4ThreeVector XY;
	double mT = (B.y()-A.y())/(B.x()-A.x());
	double cT = -mT*A.x()+A.y();
	
	XY.setX((cT-c)/(m-mT));
	XY.setY( m*XY.x() +c);
	XY.setZ(A.z());
	
	return XY;
	
}

bool uRwell_strip::pointOnsegment(G4ThreeVector X, G4ThreeVector A, G4ThreeVector B){
	
	
	if((X.x()>= fmin(A.x(), B.x())) && (X.x()<=fmax(A.x(), B.x()))&&(X.y()>=fmin(A.y(), B.y())) && (X.y()<=fmax(A.y(), B.y()))){
		
		return true;
	}else{
		return false;
	}
}

double uRwell_strip::cal_length(G4ThreeVector A, G4ThreeVector B){
	double length=0;
	length = sqrt(pow((A.x()-B.x()),2) + pow((A.y()-B.y()),2)) ;
	return length;
}

double uRwell_strip::GetBinomial(double n, double p){
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

int uRwell_strip::Number_of_strip(uRwellConstants uRwellc){
	
	int N;
	// C-------------D //
	//  -------------  //
	//   -----------   //
	//    A-------B   //
	
	/*** number of strip in AB***/
	
	int n_AB = abs(2*uRwellc.Xhalf_base/(uRwellc.stripU_pitch/sin(M_PI*uRwellc.find_stereo_angle()/180)));
	
	/** number of strip in CA **/
	double AC = sqrt((pow((uRwellc.Xhalf_base-uRwellc.Xhalf_Largebase),2) + pow((uRwellc.Yhalf+uRwellc.Yhalf),2)));
	double theta = acos(2*uRwellc.Yhalf/(AC));
	int n_AC = AC/(uRwellc.stripU_pitch/cos(theta-abs(M_PI*uRwellc.find_stereo_angle())/180));
	
	N = n_AB + n_AC+1;
	return N;
}

int uRwell_strip::strip_id(int i, uRwellConstants uRwell ){
	int ID = 0;
	// int strip=0;
	
	ID = Number_of_strip(uRwell)/2+i;
	
	return ID;
	
}


