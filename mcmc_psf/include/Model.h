#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include "Img.h"
#include "nr3.h"

class Model {

  public:
  Model (Img&, unsigned);
  long nx;
  long ny;
  long nf;

  int nedge;

  MatDoub* par_low;
  MatDoub* par_low_gs;
  MatDoub* raw_imgs;
  MatDoub* curr_imgs;
  MatDoub* gs_imgs;
  VecDoub hpar;
  VecDoub hpar_gs;

  double mp;
  double mu;
  double c;
  double freq0;
  double kapa;
  double h;
  double hkc;
  VecDoub freqr;
  double beamsize;
  double d;
  double barea;
  VecDoub p1;
  double p2; 

  VecDoub mk_raw_img (Img&, MatDoub*, int, int);
  double cal_lglike_low (Img&, int ix, int iy, int flag); //flag=0, cal likelihood of current image //flag=1, cal likelihood of guess image
  double cal_cov(Img& img, int flag);
  double cal_cov(Img& img, int ix, int iy, int flag);
  bool if_sgra(Img& img, int ix, int iy);
  bool if_bk(Img& img, int ix, int iy);
  bool if_lowb(Img& img, int ix, int iy);


}; 

#endif
