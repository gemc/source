// gemc headers
#include "ftm_strip.h"
#include "Randomize.hh"
#include <iostream>
#include <cmath>



vector<double> ftm_strip::FindStrip(int layer, double x, double y, double z, double Edep, detector Detector,ftmConstants ftmcc)
{
    // the return vector is always in pairs.
    // The first number is the ID,
    // the second number is the sharing percentage
    vector<double> strip_id;
    // number of electrons (Nt)
    int Nel = (int) (1e6*Edep/ftmcc.w_i);
    
    double r = sqrt(x*x+y*y);

    // get layer position and drift distance from the detector dimensions
    double z0           = -Detector.dimensions[2];
    double hdrift       =  Detector.dimensions[2]*2;
//    cout << Detector.name << " " << Detector.dimensions.size();
//    for(int i=0; i<Detector.dimensions.size(); i++) cout << " " << Detector.dimensions[i];
//    cout << endl;
//    cout << " z0 " <<z0 << " drift " << hdrift << " z " << z << endl;
//    cout << " layer " << layer << " x " << x << " y " << y << " z " << z << " nel " << Nel << endl;
    // old if(fabs(z-Z0[layer])>(hDrift+0.2)) cout << "Warning! z position of the FTM hit is not in the sensitive volume: " << z-Z0[layer]<< endl;
    if(fabs(z-z0)>(hdrift+0.2)) cout << "Warning! z position of the FTM hit is not in the sensitive volume: " << z-z0<< endl;
    if(r<Detector.dimensions[0] || r>Detector.dimensions[1]) cout << "Warning! r position of the FTM hit is not in the sensitive volume: " << r << endl;
  
    double sigma_td_max = ftmcc.sigma_0*sqrt(hdrift/cm);
    double sigma_td     = sigma_td_max* sqrt(fabs(z-z0)/hdrift); // expression without Lorentz angle
  
    // check if hit is within active area
    
    if(r<ftmcc.rmax && r>ftmcc.rmin) {
        double x_real=x;
        double y_real=y;
        if (layer%2==0){ // y_real is the coordinate given by the detector; strips are along x_real
            x_real = y;
            y_real = x;
        }
        if (layer%2==1){  //y_real is the coordinate given by the detector; strips are along x_real
            x_real = x;
            y_real = y;
        }
//        cout << " nel " << Nel << " x_real " << x_real << " y_real " << y_real << " z " << z << endl;
        
        
        // calculate the y coordinate of closest strip (N.B. y_strip may exceed the actual strip range
        // because the active area is a bit larger than the one instrumented with strips)
        double y_strip = ftmcc.pitch*(floor(y_real/ftmcc.pitch))+ftmcc.pitch/2;
//        cout << " nel " << Nel << " x_real " << x_real << " y_strip " << y_strip << " z " << z << endl;
        
        // calculates y coordinates of strips that may be involved in the cluster
        double weight_tot=0;
        for(int i=-ftmcc.nb_sigma; i<=ftmcc.nb_sigma; i++) {
            
            // calculate strip coordinate
            double iy = y_strip+i*ftmcc.pitch;
            
            // determine the strip number
            int istrip = get_strip_ID(x_real,iy,ftmcc);
            // check if strip number is well defined (i.e. now really inside the active volume)
            if(istrip>0) {
                // calculate x as mid point of the strip and the strip length
//                double ix = get_strip_X(x_real,iy,ftmcc);
                // calculate weight
                double weight= exp(-(iy-y_real)*(iy-y_real)/2/sigma_td/sigma_td);
                int weight_n = round(weight*Nel*10);   // round weight to closer integer to get rid of very low energy hits
                if(weight_n>0 || i==0) {               // save strips with non zero weight or at least the central one
                    strip_id.push_back(istrip);
                    strip_id.push_back(weight);
                    weight_tot += weight;
                }
            }
        }
//        cout << strip_id.size() << endl;
        if(strip_id.size()<2) {    // if not strips have been identified, set at least 1
            strip_id.push_back(-1);
            strip_id.push_back(1);
        }
        // normalize weight to have integral equal to 1
        for (unsigned i=0; i<strip_id.size()/2; i++) {
            strip_id[i*2+1]=strip_id[i*2+1]/weight_tot;
//            cout << strip_id[i*2] << " " << strip_id[i*2+1] << endl;
        }
    }
    else {  //  return strip id -1 is hit is outside active area
        strip_id.push_back(-1);
        strip_id.push_back(1);
    }
 
  return strip_id;
}


int ftm_strip::get_strip_ID(double x, double y,ftmConstants ftmcc) {
    double r=sqrt(x*x+y*y);
    if(r<ftmcc.rmax && r>ftmcc.rmin && fabs(y)<ftmcc.pitch*ftmcc.nstrips*2/6) {
        int strip = (int) floor(y/ftmcc.pitch) + 1 + ftmcc.nstrips*2/6;
        if(strip>ftmcc.nstrips*3/6) { // strip in the top sector
            strip += ftmcc.nstrips*2/6;
        }
        else if(strip>ftmcc.nstrips/6 && x>0){
            strip += ftmcc.nstrips*2/6;
        }
        return strip;
    }
    else {
        return -1;
    }
}

double ftm_strip::get_strip_X(double x, double y,ftmConstants ftmcc) {
    int istrip = get_strip_ID(x,y,ftmcc);
    double strip_X = 0;
    if(istrip<ftmcc.nstrips/6 || istrip>ftmcc.nstrips*5/6) {
        strip_X = 0;
    }
    else {
        strip_X=(ftmcc.rmax*sin(acos(fabs(y)/ftmcc.rmax))+ftmcc.rmin*sin(acos(fabs(y)/ftmcc.rmin)))/2;
        if(x<0) strip_X = -strip_X;
    }
    return strip_X;
}

/*
double ftm_strip::get_strip_L(double x, double y,ftmConstants ftmcc) {
    int istrip = get_strip_ID(x,y,ftmcc);
    double strip_L = 0;
    if(istrip<ftmcc.nstrips/6 || istrip>ftmcc.nstrips*5/6) {
        strip_L = ftmcc.rmax*sin(acos(fabs(y)/ftmcc.rmax));
    }
    else {
        strip_L = ftmcc.rmax*sin(acos(fabs(y)/ftmcc.rmax))-ftmcc.rmin*sin(acos(fabs(y)/ftmcc.rmin));
    }
    return strip_L;
}
 */

