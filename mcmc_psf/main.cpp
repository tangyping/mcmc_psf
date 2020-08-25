#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "Img.h"
#include "Model.h"
#include "Sampler.h"
#include "nr3.h"
#include "Writer.h"
#include "omp.h"

double mid(VecDoub arr, int size){

  int i,t,f;
  double even,odd;
  t=(size/2);
  f=(size/2)-1;

  even=(((arr[t]+arr[f]))/2);
  odd=arr[t];
  double val=0.;
  for(i=0;i<size;i++){
    if(size%2==0) val=even;
    if(size%2==1) val=odd;
  }
  return val;
}

int main () {

  /////////setup up begin//////
  
  int nf = 5; //number of images/bands
  
  string fsig[] = {"signal_1.fits", "signal_2.fits", "signal_3.fits", "signal_4.fits", "signal_5.fits"}; //signal maps
  vector<string> fsigs(fsig, fsig+nf);
  string ferr[] = {"err_1.fits", "err_2.fits", "err_3.fits", "err_4.fits", "err_5.fits"}; //error maps
  vector<string> ferrs(ferr, ferr+nf);
  string fker[] = {"ker_1.fits", "ker_2.fits", "ker_3.fits", "ker_4.fits", "ker_5.fits"}; //PSF kernels
  vector<string> fkers(fker, fker+nf);
  string lb[] = {"l_arr.fits", "b_arr.fits"}; //Galactic longitudes and latitudes
  vector<string> lbs(lb, lb+2);

  double wave[] = {0.016, 0.025, 0.035, 0.05, 0.11}; //wavelengths, in centimeter

  double lows[] = {20.5, log(10.), 0.5}; //lower and upper limits
  double ups[] = {24.5, log(150.), 3.5};
  VecDoub lower(3);
  VecDoub upper(3);
  for(int i=0;i<3;i++) {
    lower[i] = lows[i];
    upper[i] = ups[i];
  }

  long nburn = 200; //number of burn-in steps
  long lavg = 200; //number of sampling steps
  long nrun = nburn + lavg;

  ///parallelization control
  int nxtrd=2; //number of threads in x-axis
  int nytrd=2; //number of threads in y-axis
  /// notice that img.nx/nxtrd and img.ny/nytrd must be integers, which also need to be larger than max(img.nx_ker) and max(img.ny_ker), respectively.

  int if_hie=0;
  
  /////////end setup/////////////
  
  VecDoub waves(nf);
  for(int i=0;i<nf;i++) waves[i] = wave[i];

  Img img(fsigs, ferrs, fkers, lbs, waves);
  Model model(img, 0);
  Sampler sampler(model);

  MatDoub* raw_imgs = new MatDoub[model.nf];
  for(int i=0; i<model.nf; i++){
    raw_imgs[i] = MatDoub(model.nx, model.ny);
  }
  
  cout << img.nx << "," << img.ny << endl;

  // VecDoub as(nrun);

  VecDoub* p1_curr = new VecDoub[img.nx*img.ny];
  VecDoub* p2_curr = new VecDoub[img.nx*img.ny];
  VecDoub* p3_curr = new VecDoub[img.nx*img.ny];

  for(int i=0; i<img.nx; i++){
    for(int j=0; j<img.ny; j++){
      p1_curr[i*img.ny+j] = VecDoub(lavg);
      p2_curr[i*img.ny+j] = VecDoub(lavg);
      p3_curr[i*img.ny+j] = VecDoub(lavg);
    }
  }

  VecDoub p1b = VecDoub(img.nx*img.ny);
  VecDoub p2b = VecDoub(img.nx*img.ny);
  VecDoub p3b = VecDoub(img.nx*img.ny);

  unsigned seed=1000;

  sampler.initial_raw(img, &model);

  for(int i=0; i<nrun; i++){

    sampler.istep = sampler.istep+1;
    int ipx=20;
    int ipy=20;

    //show local parameters at each step
    cout << i << "," << model.par_low[0][ipx][ipy] << "," << 
    model.par_low[1][ipx][ipy] << "," << model.par_low[2][ipx][ipy] << endl;

    //show hyper parameters at each step
    //cout << i << "," << model.hpar[0] << "," << model.hpar[1] << "," << model.hpar[2] << "," << model.hpar[3] << "," << model.hpar[4] << endl;

    int ntrd = nxtrd*nytrd;

    int xstep=img.nx/nxtrd;
    int ystep=img.ny/nytrd;


#pragma omp parallel shared(model, seed)
{

  for(int j=0; j<img.nx/nxtrd; j++){
    for(int k=0; k<img.ny/nytrd; k++){
      #pragma omp for schedule(dynamic)
      for(int itrd=0; itrd<ntrd; itrd++){
        int ixtrd = itrd % nxtrd;
        int iytrd = (itrd - ixtrd)/nxtrd;

	if(if_hie==1){
          sampler.slice_multi(img, &model, lower, upper, j+ixtrd*xstep, k+iytrd*ystep, 1, seed+itrd);
	  seed += ntrd;
	}

	if(if_hie==0){
          sampler.slice_single(img, &model, 0, lower[0], upper[0], j+ixtrd*xstep, k+iytrd*ystep, 0.5, 1000, seed+itrd, 0);
          seed += ntrd; 
          sampler.slice_single(img, &model, 1, lower[1], upper[1], j+ixtrd*xstep, k+iytrd*ystep, 0.5, 1000, seed+itrd, 0);
          seed += ntrd;        
          sampler.slice_single(img, &model, 2, lower[2], upper[2], j+ixtrd*xstep, k+iytrd*ystep, 0.5, 1000, seed+itrd, 0);
          seed += ntrd;
	}
      }
    }
  }

}

 if(if_hie==1){
   sampler.slice_cov(img, &model, 0, log(5.), log(60.), 1000, seed);
   seed += 1;
   sampler.slice_cov(img, &model, 1, 1., 3., 1000, seed);
   seed += 1;
   sampler.slice_cov(img, &model, 2, 0.01, 0.4, 1000, seed);
   seed += 1;
   sampler.slice_cov(img, &model, 3, 0.01, 0.4, 1000, seed);
   seed += 1;
   sampler.slice_cov(img, &model, 4, -1., 1., 1000, seed);
   seed += 1;
 }
 
   if(i>=nburn){ 
     for(int ix=0; ix < img.nx; ix++){
       for(int iy=0; iy < img.ny; iy++){
        p1_curr[ix*img.ny+iy][(i-nburn)%lavg] = model.par_low[0][ix][iy];
        p2_curr[ix*img.ny+iy][(i-nburn)%lavg] = model.par_low[1][ix][iy];
        p3_curr[ix*img.ny+iy][(i-nburn)%lavg] = model.par_low[2][ix][iy];
       }
     }	

   }


   if((i-nburn) % lavg == (lavg-1) && i>nburn){

     for(int ix=0; ix<img.nx; ix++){
       for(int iy=0; iy<img.ny; iy++){
         p1b[ix*img.ny+iy] = 0.;
         p2b[ix*img.ny+iy] = 0.;
         p3b[ix*img.ny+iy] = 0.;
       }
     }
		 
     for(int ix=0; ix < img.nx; ix++){
       for(int iy=0; iy < img.ny; iy++){
         p1b[ix*img.ny+iy] = mid(p1_curr[ix*img.ny+iy], lavg);
         p2b[ix*img.ny+iy] = mid(p2_curr[ix*img.ny+iy], lavg);
         p3b[ix*img.ny+iy] = mid(p3_curr[ix*img.ny+iy], lavg);
       }		
     }
		
     writeVecOut("p1_med.dat", p1b);
     writeVecOut("p2_med.dat", p2b);
     writeVecOut("p3_med.dat", p3b);

    }

  }

  unsigned long rows(3);     
  string hduName("TABLE_BINARY");          
  vector<string> colName(3,"");
  vector<string> colForm(3,"");
  vector<string> colUnit(3,"");

  colName[0] = "p1";
  colName[1] = "p2";
  colName[2] = "p3";

  ostringstream oss;
  oss << lavg;
  string str1 = oss.str();
  string str2 = "D";

  colForm[0] = str1+str2;
  colForm[1] = str1+str2;
  colForm[2] = str1+str2;

  colUnit[0] = "cm2";
  colUnit[1] = "K";
  colUnit[2] = "";  

  unique_ptr<FITS> pFits(new FITS("table.fits",Write));

  vector<valarray<double> > p1_samp(img.nx*img.ny);
  vector<valarray<double> > p2_samp(img.nx*img.ny);
  vector<valarray<double> > p3_samp(img.nx*img.ny);

  for(int ix=0; ix<img.nx; ix++){
    for(int iy=0; iy<img.ny; iy++){
      p1_samp[ix*img.ny+iy].resize(lavg);
      p2_samp[ix*img.ny+iy].resize(lavg);
      p3_samp[ix*img.ny+iy].resize(lavg);
      for(long i=0; i<lavg; i++){
	 p1_samp[ix*img.ny+iy][i] = p1_curr[ix*img.ny+iy][i];
	 p2_samp[ix*img.ny+iy][i] = p2_curr[ix*img.ny+iy][i];
	 p3_samp[ix*img.ny+iy][i] = p3_curr[ix*img.ny+iy][i];
      }
    }
  }

  Table* newTable = pFits->addTable(hduName,rows,colName,colForm,colUnit);
  newTable->column(colName[0]).writeArrays(p1_samp,0);  
  newTable->column(colName[1]).writeArrays(p2_samp,0);  
  newTable->column(colName[2]).writeArrays(p3_samp,0);  

  return 0;

}
