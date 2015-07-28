#include "bst_strip.h"
#include "fst_strip.h"

#include <iostream>
using namespace std;

int main(int argn, char** argv)
{
 if(argn != 7)
 {
    cout << endl << " Wrong mumber of arguments. Usage: example layer sector x y z BST" << endl << endl;
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

 if(D == "BST")
 {
    class bst_strip bsts;
    bsts.fill_infos();
    bsts.FindCard(layer, z);
    X = bsts.x;
    Y = bsts.y;
    Z = bsts.z;
    Strip = bsts.FindStrip(layer, sector, X, Y, Z);
 }
 if(D == "FST")
 {
    class fst_strip fsts;
    fsts.fill_infos();
    X = x;
    Y = y;
    Z = z;
    Strip = fsts.FindStrip(layer, sector, X, Y, Z);
 }


 cout << D << ": "; 
 cout << " x = " << X << " " ;
 cout << " y = " << Y << " " ;
 cout << " z = " << Z << " " ;
 cout << " Strip = " << Strip << endl;
 
 return 1;

}
