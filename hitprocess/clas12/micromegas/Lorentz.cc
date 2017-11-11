#include "Lorentz.h"
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>

#include <cmath>
#include <iostream>
using namespace ccdb;

Lorentz::Lorentz(){
  
}

Lorentz::~Lorentz(){

}

void Lorentz::Initialize(int runno){

  if(getenv ("CCDB_CONNECTION") != NULL) {
    connection = (string) getenv("CCDB_CONNECTION");
  } else {
    connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
  }

  variation  = "default";
  vector<vector<double> > data;
  auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(connection));
  
  sprintf(database,"/calibration/mvt/lorentz");
  data.clear(); calib->GetCalib(data,database);
  
  float pe=0;
  float pb=0;
   int i = 0;
  for(unsigned row = 0; row < data.size(); row++)
    {
       Lor_grid.push_back(data[row][2]);
       E_grid.push_back(data[row][0]);
       B_grid.push_back(data[row][1]);
       if( Ne==0 && Nb==0 ){
	 emin = data[row][0];
	 emax = data[row][0];
	 bmin = data[row][1];
	 bmax = data[row][1];

	 i++;
	 Ne++;
	 Nb++;
	 pe = data[row][0];
	 pb = data[row][1];
	 continue;
       }
       
       // check max and minima
       if ( data[row][0] < emin ) emin = data[row][0];
       if ( data[row][0] > emax ) emax = data[row][0];
       if ( data[row][1] < bmin ) bmin = data[row][1];
       if ( data[row][1] > bmax ) bmax = data[row][1];
       
       // count E and B
       if( fabs( pe - data[row][0] ) > 0.0001 ) Ne++;
       if ( Ne==1) { if( fabs( pb - data[row][1] ) > 0.0001 ) Nb++; } // only for the first Nb value
       
      
       pe = data[row][0];
       pb = data[row][1];
       
    }
}

int Lorentz::getBin( float e, float b){
  if( Lor_grid.size()==0 ) return 0;
  float de = (emax-emin)/(Ne-1);
  float db = (bmax-bmin)/(Nb-1);
  
  int ie = floor( (e - emin)/de );
  int ib = floor( (b - bmin)/db );
  
//   std::cout << ie << "  " << ib << "\n";
  
  return ib + Nb * ie ;
}

float Lorentz::GetAngle(float xe, float xb){
  if(Lor_grid.size()==0||xe==0||xb==0) return 0.;
  float de = (emax-emin)/(Ne-1);
  float db = (bmax-bmin)/(Nb-1);

  if (xe<emin) {
    xe=emin;
    cout<<"Warning: E out of grid... setting it to Emin"<<endl;
  }
  if (xe>=emax) {
    xe=emax*0.99;
    cout<<"Warning: E out of grid... setting it to Emax"<<endl;
  }
  if (xb>bmax) {
    xb=bmax*0.99;
    cout<<"Warning: B out of grid... setting it to Bmax"<<endl;
  }
  
  int i11 = getBin( xe, xb);
  int i12 = getBin( xe, xb+db);
  int i21 = getBin( xe+de, xb);
  int i22 = getBin( xe+de, xb+db);

  float Q11 = 0; float Q12 = 0; float Q21 = 0;   float Q22 = 0;
  float e1 = emin; float e2 = emax; float b1 = 0; float b2 = bmax; 
  if (i11>=0) {
    Q11=Lor_grid.at(i11); e1 = E_grid.at(i11);  b1 = B_grid.at(i11);
  }
  if (i12>=0) Q12 = Lor_grid.at(i12);
  if (xb>=bmin) Q21 = Lor_grid.at(i21);
  if (xb<bmin) Q21 = 0;
  if (i22>=0) {
    Q22 = Lor_grid.at(i22); e2 = E_grid.at(i22);  b2 = B_grid.at(i22);
  }
 
  float R1 = linInterp( xe, e1,e2,Q11,Q21);
  float R2 = linInterp( xe, e1,e2,Q12,Q22);
  
  float P =  linInterp( xb, b1,b2,R1, R2);
  
  return P;
}

float Lorentz::linInterp( float x, float x1, float x2, float y1, float y2 ){
  // linear interpolation
  // return y = f(x), given x1, y1=f(x1) and x2, y2=f(x2) 
  // y = m * ( x - x1 ) + y1
  // m = ( y2 - y1)/(x2 - x1)
  
  // compute m
  float m = (y2 - y1)/(x2 - x1);
  
  // return
  return m * ( x - x1 ) + y1;
}
