#define _USE_MATH_DEFINES

#include "GslRandom.h"
#include <iostream>
#include "Model.h"
#include "Img.h"
#include "Sampler.h"
#include <math.h>

using namespace std;

Sampler::Sampler(Model& model){

  istep=0; //current step
  
}

void Sampler::initial_raw(Img& img, Model* pm){

  VecDoub flux;

  for(int ix=0; ix < pm->nx; ix++){
    for(int iy=0; iy < pm->ny; iy++){

      flux = pm->mk_raw_img(img, pm->par_low, ix, iy);

      for(int iwav=0; iwav<img.nf; iwav++){
			  pm->raw_imgs[iwav][ix][iy] = flux[iwav];
			  pm->curr_imgs[iwav][ix][iy] = 0.; //current image
			  pm->gs_imgs[iwav][ix][iy] = 0.;//guess image
      }
    }
  }

  //convolution
  
  for(int iwav=0; iwav<img.nf; iwav++){

    for(long ix=0; ix<img.nx; ix++){
      for(long iy=0; iy<img.ny; iy++){

        for(long j=0; j<img.nx_ker[iwav]; j++){
          for(long k=0; k<img.ny_ker[iwav]; k++){

            int ixp = ix+j-(img.nx_ker[iwav]-1)/2;
            int iyp = iy+k-(img.ny_ker[iwav]-1)/2;
            int ixr = ixp;
            int iyr = iyp;
            if(ixp < 0) ixr = img.nx+ixp;
            if(ixp >= img.nx) ixr = ixp-img.nx;
            if(iyp < 0) iyr = img.ny+iyp;
            if(iyp >= img.ny) iyr = iyp-img.ny;
            pm->curr_imgs[iwav][ixr][iyr] = pm->curr_imgs[iwav][ixr][iyr]+img.kers[iwav][j][k] * pm->raw_imgs[iwav][ix][iy];
            pm->gs_imgs[iwav][ixr][iyr] = pm->gs_imgs[iwav][ixr][iyr]+img.kers[iwav][j][k] * pm->raw_imgs[iwav][ix][iy];
	  }
        }
      }
    }
  }

}

void Sampler::slice_multi(Img& img, Model* pm, VecDoub lower, VecDoub upper, int ix, int iy, int mtype, unsigned seed){

  GslRandom ran(seed);
  int np=3;

  VecDoub x0(np);
  VecDoub Lhat(np);
  VecDoub Rhat(np);

  double y;
  double e;

  for(int ip=0; ip<np; ip++) x0[ip] = pm->par_low[ip][ix][iy];  
  e = ran.expDeviate(1.);
  y = get_post(img, pm, ix, iy, e, mtype, 0);

  //start shrink
  
  for(int ip=0; ip<np; ip++){
    Lhat[ip] = lower[ip];
    Rhat[ip] = upper[ip];
  }

  VecDoub U(np);
  VecDoub x1(np);
  double lk_x1;

  for(int ip=0; ip<np; ip++){
    U[ip] = ran.uniformDeviate(0., 1.);
    x1[ip] = Lhat[ip] + U[ip]*(Rhat[ip] - Lhat[ip]);
    pm->par_low_gs[ip][ix][iy] = x1[ip];
  }

  lk_x1 = get_post(img, pm, ix, iy, e, mtype, 1);
  
  for(int ip=0; ip<np; ip++){
    if(x1[ip] < x0[ip]) {
      Lhat[ip] = x1[ip];
    }
    else{
      Rhat[ip] = x1[ip];
    }
  }
  
  while(y > lk_x1 && (Rhat[0] - Lhat[0] > 1e-10 || Rhat[1] - Lhat[1] > 1e-10 || Rhat[2] - Lhat[2] > 1e-10)){

    for(int ip=0; ip<np; ip++){
      if(Rhat[ip] - Lhat[ip] > 1e-10){
        U[ip] = ran.uniformDeviate(0., 1.);
        x1[ip] = Lhat[ip] + U[ip]*(Rhat[ip] - Lhat[ip]);
        pm->par_low_gs[ip][ix][iy] = x1[ip];
      }
    }

    lk_x1 = get_post(img, pm, ix, iy, e, mtype, 1);

    for(int ip=0; ip<np; ip++){
      if(Rhat[ip] - Lhat[ip] > 1e-10){
        if(x1[ip] < x0[ip]) {
          Lhat[ip] = x1[ip];
        }
        else{
          Rhat[ip] = x1[ip];
	}
      }
    }
  }

  //update parameter
  
  for(int ip=0; ip<np; ip++) pm->par_low[ip][ix][iy] = x1[ip];

  //update current image
  
  for(int iwav=0; iwav<img.nf; iwav++){
    for(long j=0; j<img.nx_ker[iwav]; j++){
      for(long k=0; k<img.ny_ker[iwav]; k++){
        int ixp = ix+j-(img.nx_ker[iwav]-1)/2;
        int iyp = iy+k-(img.ny_ker[iwav]-1)/2;
        int ixr = ixp;
        int iyr = iyp;
        if(ixp < 0) ixr = img.nx+ixp;
        if(ixp >= img.nx) ixr = ixp-img.nx;
        if(iyp < 0) iyr = img.ny+iyp;
        if(iyp >= img.ny) iyr = iyp-img.ny;        
        pm->curr_imgs[iwav][ixr][iyr] = pm->gs_imgs[iwav][ixr][iyr];
      }
    }
  }

}


void Sampler::slice_single(Img& img, Model* pm, int ip, double lower, double upper, int ix, int iy, double w, int m, unsigned seed, int mtype){

  GslRandom ran(seed);

  double x0;
  double Lhat;
  double Rhat;
  double e;
  double y;

  x0 = pm->par_low[ip][ix][iy];
  e = ran.expDeviate(1.);
  y = get_post(img, pm, ix, iy, e, mtype, 0);

  //stepping out
  /*
  double J;
	double K;
	double V;
	double fL;
	double fR;
  double V;
	
  L = pm->par_low[ip][ix][iy] - w*u;
	R = L+w;

	V = ran.uniformDeviate(0.,1.);
	J = floor(m*V);
	K = (m-1) - J;

	pm->par_low_gs[ip][ix][iy] = L;
	fL = get_post(img, pm, ix, iy, e, mtype, 1);

	pm->par_low_gs[ip][ix][iy] = R;
	fR = get_post(img, pm, ix, iy, e, mtype, 1);

  while(J > 0 && y < fL && L > lower){
		L=L-w;
		J=J-1;
		pm->par_low_gs[ip][ix][iy]=L;
		fL = get_post(img, pm, ix, iy, e, mtype, 1);
	}
	
  while(K > 0 && y < fR && R < upper){
		R=R+w;
		K=K-1;
		pm->par_low_gs[ip][ix][iy]=R;
		fR = get_post(img, pm, ix, iy, e, mtype, 1);
	}

  //now start shrinkage
  Lhat = L;
  Rhat = R;

  if(Lhat<lower) Lhat=lower;
  if(Rhat>upper) Rhat=upper;
  */

  //shrink 
  Lhat=lower; //directly shrink from upper & lower limits (no stepping out)
  Rhat=upper;
	
  double U;
  U = ran.uniformDeviate(0., 1.);
  double x1;
  x1 = Lhat + U*(Rhat - Lhat);
  double lk_x1;
  pm->par_low_gs[ip][ix][iy] = x1;

  lk_x1 = get_post(img, pm, ix, iy, e, mtype, 1);

  if(x1 < x0) {
    Lhat = x1;
  }
  else{
    Rhat = x1;
  }
  
  while(y > lk_x1 && Rhat - Lhat > 1e-10){

    U = ran.uniformDeviate(0., 1.);
    x1 = Lhat + U*(Rhat - Lhat);
    pm->par_low_gs[ip][ix][iy] = x1;

    lk_x1 =  get_post(img, pm, ix, iy, e, mtype, 1);

    if(x1 < x0) {
      Lhat = x1;
    }
    else{
      Rhat = x1;
    }
  }

  //update parameter
  pm->par_low[ip][ix][iy] = x1;

  //update current image
  for(int iwav=0; iwav<img.nf; iwav++){
    for(long j=0; j<img.nx_ker[iwav]; j++){
      for(long k=0; k<img.ny_ker[iwav]; k++){
        int ixp = ix+j-(img.nx_ker[iwav]-1)/2;
        int iyp = iy+k-(img.ny_ker[iwav]-1)/2;
        int ixr = ixp;
        int iyr = iyp;
        if(ixp < 0) ixr = img.nx+ixp;
        if(ixp >= img.nx) ixr = ixp-img.nx;
        if(iyp < 0) iyr = img.ny+iyp;
        if(iyp >= img.ny) iyr = iyp-img.ny;
        pm->curr_imgs[iwav][ixr][iyr] = pm->gs_imgs[iwav][ixr][iyr];
      }
    }
  }
  
}

////////Sampler of the hyper parameters
void Sampler::slice_cov(Img& img, Model* pm, int ip, double lower, double upper, int m, unsigned seed){ //w: step size; m: upper limit of slice length

  GslRandom ran(seed);

  double x0;
  double Lhat;
  double Rhat;
  double y;
  double e;

  x0 = pm->hpar[ip];
  e = ran.expDeviate(1.);
  y = pm->cal_cov(img, 0) - e;

  Lhat = lower;
  Rhat = upper;
    
  double U;
  U = ran.uniformDeviate(0., 1.);
  double x1;
  x1 = Lhat + U*(Rhat - Lhat);
  double lk_x1;
  pm->hpar_gs[ip] = x1;
  lk_x1 = pm->cal_cov(img, 1);
  if(x1 < x0){
    Lhat = x1;
  }
  else{
    Rhat = x1;
  }

  while(y > lk_x1 && Rhat - Lhat > 1e-10){

    U = ran.uniformDeviate(0., 1.);
    x1 = Lhat + U*(Rhat - Lhat);
    pm->hpar_gs[ip] = x1;
    lk_x1 = pm->cal_cov(img, 1);
    if(x1 < x0) {
      Lhat = x1;
    }
    else{
      Rhat = x1;
    }
  }

  pm->hpar[ip] = x1;
 
}

////get posterior////
double Sampler::get_post(Img& img, Model* pm, int ix, int iy, double e, 
                           int mtype, int flag){

    double val;
    if(flag==1) e=0.;

    if(mtype == 0){
      val = pm->cal_lglike_low(img, ix, iy, flag) - e;
    }
    if(mtype == 1){
      val = pm->cal_lglike_low(img, ix, iy, flag) + 
	    pm->cal_cov(img, ix, iy, flag) - e;
    }

    return val;

}

////////////
