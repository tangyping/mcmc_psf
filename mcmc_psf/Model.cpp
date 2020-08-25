#define _USE_MATH_DEFINES
#include "GslRandom.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include "Model.h"
#include "Img.h"
#include <math.h>
#include "omp.h"

using namespace std;

//// The constructor just initializes parameters //
Model::Model(Img& img, unsigned seed){

  nx = img.nx;
  ny = img.ny;
  nf = img.nf;

  //unconvolved model-images
  raw_imgs = new MatDoub[img.nf];
  for(int i=0; i<img.nf; i++){
    raw_imgs[i] = MatDoub(nx, ny);
  }

  //current model-images
  curr_imgs = new MatDoub[img.nf];
  for(int i=0; i<img.nf; i++){
    curr_imgs[i] = MatDoub(nx, ny);
  }

  //guess model-images
  gs_imgs = new MatDoub[img.nf];
  for(int i=0; i<img.nf; i++){
    gs_imgs[i] = MatDoub(nx, ny);
  }

  int npar = 3;

  par_low = new MatDoub[npar];
  for(int i=0; i<npar; i++){
    par_low[i] = MatDoub(img.nx, img.ny);
  }

  par_low_gs = new MatDoub[npar];
  for(int i=0; i<npar; i++){
    par_low_gs[i] = MatDoub(nx, ny);
  }

  ////set up random generator////
  GslRandom ran(seed);

  ////Now initalizing parameters ////
  for(int i=0; i<img.nx; i++){
    for(int j=0; j<img.ny; j++){

      par_low[0][i][j] =  ran.uniformDeviate(22.0, 23.0); // Column Density
      par_low_gs[0][i][j] = par_low[0][i][j];
      par_low[1][i][j] =  ran.uniformDeviate(log(20.), log(30.)); // Temperature 30
      par_low_gs[1][i][j] = par_low[1][i][j];
      par_low[2][i][j] =  ran.uniformDeviate(1.5, 2.5); // Beta
      par_low_gs[2][i][j] = par_low[2][i][j];

    }
  }

  for(int i=0; i<img.nf; i++){
    for(int j=0; j<img.nx; j++){
      for(int k=0; k<img.ny; k++){
        curr_imgs[i][j][k] = 0.;
        gs_imgs[i][j][k] = 0.;
      }
    }
  }

  //initial guess of hyper parameters (convariance matrix)
  Doub hm_i[] = { log(20.), 1.8 }; //means 
  Doub hs_i[] = { 0.2, 0.2}; //sigmas

  hpar = VecDoub(5);
  hpar_gs = VecDoub(5);

  for(int i=0; i<2; i++) {
    hpar[i] = hm_i[i]; //means
    hpar_gs[i] = hm_i[i];
    hpar[i+2] = hs_i[i]; //sigmas
    hpar_gs[i+2] = hs_i[i];
  }

  hpar[4] = 0.; //correlation cofficient
  hpar_gs[4] = 0.;

  
  //Planck function
  mp = 1.6726e-24;
  mu = 2.8;
  c = 2.9979e10; //cm
  freq0 = c/0.1;
  kapa = 1.37;
  h = 6.62608e-27;
  hkc = 1.43875;
  freqr = VecDoub(img.nf);
  for(int i=0; i<img.nf; i++){
    freqr[i] = c/img.wave[i]/freq0;
  } 

  beamsize = 9.5/2.355;
  barea = 2. * M_PI * pow(beamsize/3600. * M_PI / 180., 2.);
	
  p1 = VecDoub(img.nf);
  for(int i=0; i<img.nf; i++){
    p1[i] = 2.* h * pow(c/img.wave[i], 3);
  }
  p2 = pow(c,2);

}


VecDoub Model::mk_raw_img(Img& img, MatDoub* par, int ix, int iy){  
//The dimension of par should match img
//only change the value of pixel(ix,iy) in curr_imgs

  VecDoub val(img.nf);

  double Ncol;
  Ncol = pow(10., par[0][ix][iy]);

  for(int i=0; i<img.nf; i++){
        val[i] = p1[i]/(p2 * ( exp( hkc/(exp(par[1][ix][iy]) *img.wave[i])) -1. )) * (1. - exp(-kapa * pow(freqr[i], par[2][ix][iy]) * Ncol * mp * mu * 0.01)) * barea * 1e23;
  }

  return val;
  
}

//calculate likelihood/the change of likelihood due to updating of parameters at (ix,iy)
double Model::cal_lglike_low(Img& img, int ix, int iy, int flag){
 
  double lglike=0.;
  VecDoub flux1;
  VecDoub flux2;

  int ixp=0;
  int iyp=0;
  int ixr=0;
  int iyr=0;
  
  if(flag==0) {
    for(int iwav=0; iwav<img.nf; iwav++){
      for(long j=0; j<img.nx_ker[iwav]; j++){
        for(long k=0; k<img.ny_ker[iwav]; k++){
          ixp = ix+j-(img.nx_ker[iwav]-1)/2;
          iyp = iy+k-(img.ny_ker[iwav]-1)/2;
          ixr = ixp;
          iyr = iyp;
          if(ixp <0 or ixp >=img.nx) {
            if(ixp < 0) ixr = img.nx+ixp;
            if(ixp >= img.nx) ixr = ixp-img.nx;
	  }
          if(iyp <0 or iyp >=img.ny) {
            if(iyp < 0) iyr = img.ny+iyp;
            if(iyp >= img.ny) iyr = iyp-img.ny;
	  }
	  
	  if(img.signals[iwav][ixr][iyr] >= 3.*img.errs[iwav][ixr][iyr]){
	     lglike = lglike-pow(img.signals[iwav][ixr][iyr] - curr_imgs[iwav][ixr][iyr], 2)/(2.*pow(img.errs[iwav][ixr][iyr],2));	        
	  }
          else{
  	     if(curr_imgs[iwav][ixr][iyr] > 3.*img.errs[iwav][ixr][iyr]) lglike -= 1e20;
             if(curr_imgs[iwav][ixr][iyr] < 3.*img.errs[iwav][ixr][iyr]) lglike = lglike;
	  }
	}
      }
    }
  }
  
  if(flag==1){
    flux1 = mk_raw_img(img, par_low_gs, ix, iy);
    flux2 = mk_raw_img(img, par_low, ix, iy);
    for(int iwav=0; iwav<img.nf; iwav++){

      for(long j=0; j<img.nx_ker[iwav]; j++){
        for(long k=0; k<img.ny_ker[iwav]; k++){
          ixp = ix+j-(img.nx_ker[iwav]-1)/2;
          iyp = iy+k-(img.ny_ker[iwav]-1)/2;
          ixr = ixp;
          iyr = iyp;
          if(ixp <0 or ixp >=img.nx) {
            if(ixp < 0) 
              ixr = img.nx+ixp;
            if(ixp >= img.nx) 
              ixr = ixp-img.nx;
	  }
          if(iyp <0 or iyp >=img.ny) {
            if(iyp < 0) 
              iyr = img.ny+iyp;
            if(iyp >= img.ny) 
              iyr = iyp-img.ny;
	  }
          gs_imgs[iwav][ixr][iyr] = curr_imgs[iwav][ixr][iyr]+img.kers[iwav][j][k]*(flux1[iwav]-flux2[iwav]);
          if(img.signals[iwav][ixr][iyr] >= 3.*img.errs[iwav][ixr][iyr]){
	    lglike = lglike-pow(img.signals[iwav][ixr][iyr] - gs_imgs[iwav][ixr][iyr], 2)/(2.*pow(img.errs[iwav][ixr][iyr],2));
	  }					
          else{
             if(gs_imgs[iwav][ixr][iyr] > 3.*img.errs[iwav][ixr][iyr]) lglike -= 1e20;
             if(gs_imgs[iwav][ixr][iyr] < 3.*img.errs[iwav][ixr][iyr]) lglike = lglike;
	  }    
	}
      }
    }
  }

  ////smoothness prior
  double pcen;
  int nsp=3;
  VecInt isp(nsp);
  VecDoub ssig(nsp);
  isp[0] = 0;
  isp[1] = 1;
  isp[2] = 2;
  ssig[0] = 0.1;
  ssig[1] = 0.;
  ssig[2] = 0.2; //50.; //// 1/(2.sscal) = sig(dp)^2, sscal = 1/(2.*sig(dp)^2)

  if(ix>=0 && (ix<(img.nx-0)) && iy>=0 && (iy<(img.ny-0))){  
    for(int i=0; i<nsp; i++){
	
      int iip = isp[i];

      if(flag==0) pcen = par_low[iip][ix][iy];
      if(flag==1) pcen = par_low_gs[iip][ix][iy];
      int il=ix-1;
      if(il<0) il=0;
      int ir=ix+1;
      if(ir>=img.nx) ir=img.nx-1;
      int id=iy-1;
      if(id<0) id=0;
      int iu=iy+1;
      if(iu>=img.ny) iu=img.ny-1;
		
      double vl = pcen - par_low[iip][il][iy];
      double vr = pcen - par_low[iip][ir][iy];
      double vu = pcen - par_low[iip][ix][iu];
      double vd = pcen - par_low[iip][ix][id];

      if(ssig[i]>0.) lglike = lglike - 1./2.*(pow(vl,2) + pow(vr,2) + pow(vu,2) + pow(vd,2))/(2.*pow(ssig[i],2));

    }
  }

  return lglike;
}

/////////calculate the hyper-prior from the entire grids 
double Model::cal_cov(Img& img, int flag){

  double val=0.;
  double deter = pow(hpar[2]*hpar[3], 2) * (1.-pow(hpar[4],2));
  double deter_gs = pow(hpar_gs[2]*hpar_gs[3], 2) * (1.-pow(hpar_gs[4],2));
  double A;

  if(flag==0){
    for(int ix=0; ix<img.nx; ix++){
      for(int iy=0; iy<img.ny; iy++){
        if(if_sgra(img, ix, iy)==false && if_lowb(img,ix,iy)==true){//excluding sgr* and high latitudes regions

          A = pow(hpar[2],2) * pow(par_low[2][ix][iy] - hpar[1], 2) +
              pow(hpar[3],2) * pow(par_low[1][ix][iy] - hpar[0], 2) -
	      2.*hpar[4]*hpar[2]*hpar[3]*(par_low[1][ix][iy]-hpar[0])*(par_low[2][ix][iy]-hpar[1]);
        
          val += -0.5*log(deter) - 5.*log((1. + 1./(8.*deter) * A));
					
         }
       }
     }
   }

  if(flag==1){
    for(int ix=0; ix<img.nx; ix++){
      for(int iy=0; iy<img.ny; iy++){
        if(if_sgra(img, ix, iy)==false && if_lowb(img,ix,iy)==true){

          A =  pow(hpar_gs[2],2) * pow(par_low[2][ix][iy] - hpar_gs[1], 2) +
               pow(hpar_gs[3],2) * pow(par_low[1][ix][iy] - hpar_gs[0], 2) -
	       2.*hpar_gs[4]*hpar_gs[2]*hpar_gs[3]*(par_low[1][ix][iy]-hpar_gs[0])*(par_low[2][ix][iy]-hpar_gs[1]);
        
          val += -0.5*log(deter_gs) - 5.*log((1. + 1./(8.*deter_gs) * A));

	 }
       }
     }
   }

   val -= pow(log(hpar[2]) - log(3.2), 2.)/1e4  + pow(log(hpar[3]) - log(2.0), 2.)/1e4;
	
   return val;

}


/////////calculate the hyper-prior/change of hyper-prior due to updating at (ix,iy) 
double Model::cal_cov(Img& img, int ix, int iy, int flag){

  double val=0.;
  double deter = pow(hpar[2]*hpar[3], 2) * (1.-pow(hpar[4],2));
  double A;

  if(flag==0){
    if(if_sgra(img, ix, iy)==false && if_lowb(img,ix,iy)==true){

       A = pow(hpar[2],2) * pow(par_low[2][ix][iy] - hpar[1], 2) +
           pow(hpar[3],2) * pow(par_low[1][ix][iy] - hpar[0], 2) -
	   2.*hpar[4]*hpar[2]*hpar[3]*(par_low[1][ix][iy]-hpar[0])*(par_low[2][ix][iy]-hpar[1]);
        
       val = -5.*log((1. + 1./(8.*deter) * A));
			 
    }
  }

  if(flag==1){
    if(if_sgra(img, ix, iy)==false && if_lowb(img,ix,iy)==true){

       A = pow(hpar[2],2) * pow(par_low_gs[2][ix][iy] - hpar[1], 2) +
           pow(hpar[3],2) * pow(par_low_gs[1][ix][iy] - hpar[0], 2) -
	   2.*hpar[4]*hpar[2]*hpar[3]*(par_low_gs[1][ix][iy]-hpar[0])*(par_low_gs[2][ix][iy]-hpar[1]);

       val = -5.*log((1. + 1./(8.*deter) * A));        
     }
  }

  return val;

}

/////

bool Model::if_sgra(Img& img, int ix, int iy){ 

  bool val;
  val = false;
  if( abs(img.lbs[0][ix][iy] - 359.9442) < 0.0167 && 
      abs(img.lbs[1][ix][iy] - (-0.0462)) < 0.0167 ){
    val = true;
    return val;
  }

  return val;
}

///////////

bool Model::if_bk(Img& img, int ix, int iy){ 

  bool val;
  val = false;
  if( (img.lbs[0][ix][iy] > 0.15) && (img.lbs[0][ix][iy] < 0.3) &&
			(img.lbs[1][ix][iy] > 0.24) ){
    val = true;
    return val;
	}
  if( (img.lbs[0][ix][iy] >  359.58) && (img.lbs[0][ix][iy] <  359.81) &&
			(img.lbs[1][ix][iy] > 0.24) ){
    val = true;
    return val;
	}
  if( (img.lbs[0][ix][iy] >  359.45) && (img.lbs[0][ix][iy] <  359.84) &&
			(img.lbs[1][ix][iy] < -0.42) ){
    val = true;
    return val;
  }

  return val;

}

///////////////////////
bool Model::if_lowb(Img& img, int ix, int iy){ 

  bool val;
  val = false;
  if( img.lbs[1][ix][iy] < 0.1 && img.lbs[1][ix][iy] > -0.2){
    val = true;
    return val;
  }

  return val;

}
