#ifndef IMG_H
#define IMG_H

#include <CCfits/CCfits>
#include <iostream>
using namespace CCfits;
using namespace std;
//#include <vector>
#include <string>
#include "nr3.h"

class Img {    
  public:
  Img (vector<string>&, vector<string>&, vector<string>&, vector<string>&, VecDoub&);

  // file names of signal .fits
  vector<string> fn1;
  // file names of err .fits
  vector<string> fn2;
  // file names of ker .fits
  vector<string> fn3;
  // file names of tabulated sed
  vector<string> fn4; // lb array

  //number of files
  int nf;
  //number of pars
  int np;
  //wavelengths (in cm)
  VecDoub wave;
  int nx;
  int ny;

  VecInt nx_ker;
  VecInt ny_ker;
  VecInt ndim;

  //array of signal images
  MatDoub* signals;
  //array of err images
  MatDoub* errs;
  //array of kernel images;
  MatDoub* kers;
  VecDoub* tab_par;
  VecDoub* tseds;
  MatDoub* lbs;

};

#endif
