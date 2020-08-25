#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include "Img.h"
#include "Model.h"
#include "nr3.h"

class Sampler {

  public:
  Sampler (Model&);
  Model* pmodel;

  long istep;
  
  void initial_raw(Img& img, Model* pm);
  void slice_multi(Img& img, Model* pm, VecDoub lower, VecDoub upper, int ix, int iy, int m, unsigned seed);
  void slice_single(Img& img, Model* pm, int ip, double lower, double upper, int ix, int iy, double w, int m, unsigned seed, int type);
  void slice_cov(Img& img, Model* pm, int ip, double lower, double upper, int m, unsigned seed);
  double get_post(Img& img, Model* pm, int ix, int iy, double e, int mtype, int flag);
 
}; 


#endif
