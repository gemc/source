#include "fmt_strip.h"
#include "bmt_strip.h"

#include <iostream>
using namespace std;

int main(int argn, char** argv)
{
 if(argn != 7)
 {
    cout << endl << " Wrong mumber of arguments. Usage: example layer sector x y z BMT" << endl << endl;
    exit(0);
 }

 int layer  = atoi(argv[1]);
 int sector = atoi(argv[2]);
 double x   = atof(argv[3]);
 double y   = atof(argv[4]);
 double z   = atof(argv[5]);
 string D = argv[6];

 double X,Y,Z;
 
 int Strip;

 if(D == "FMT")
 {
    class fmt_strip fmts;
    fmts.fill_infos();
    X = x;
    Y = y;
    Z = z;
    Strip = fmts.FindStrip(layer, sector, X, Y, Z);
 }

 if(D == "BMT")
 {
    class bmt_strip bmts;
    bmts.fill_infos();
    X = x;
    Y = y;
    Z = z;
    Strip = bmts.FindStrip(layer, sector, X, Y, Z);
 }


 cout << D << ": "; 
 cout << " x = " << X << " " ;
 cout << " y = " << Y << " " ;
 cout << " z = " << Z << " " ;
 cout << " Strip = " << Strip << endl;
 
 return 1;

}
