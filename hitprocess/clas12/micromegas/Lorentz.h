#ifndef LORENTZ_H
#define LORENTZ_H 1

// c++
#include <vector>
#include "HitProcess.h"
#include <string>
using namespace std;


class Lorentz
{
 public:
  Lorentz(); 
  ~Lorentz();
  vector<double> Lor_grid;
  vector<double> E_grid;
  vector<double> B_grid;

  float emin = 999999.;
  float emax = 0.;
  float bmin = 999999.;
  float bmax = 0;
  int Ne = 0;
  int Nb = 0;
  
  string variation;
  string date;
  string connection;
  char   database[80];
  
  void Initialize(int runno);
  int getBin( float e, float b);
  float GetAngle(float xe, float xb);
  float linInterp( float x, float x1, float x2, float y1, float y2 );
};

#endif
